# The much of the code in the sequencing part is adapted from https://github.com/caelan/pb-construction

import Rhino.Geometry.Point3d as Point3d
import Grasshopper.DataTree as DataTree # Add, AddRange(list, GHPath(i, j))
import Grasshopper.Kernel.Data.GH_Path as GHPath
import System

from ghpythonlib.treehelpers import list_to_tree, tree_to_list

##############################

import heapq
import random
import time
import math
import copy
from collections import namedtuple, OrderedDict, defaultdict, deque

INF = 1e12
DEFAULT_TRANS_TOL = 5e-2 # meter

##############################
import os
import clr
# https://ironpython.net/documentation/dotnet/dotnet.html

clr.AddReferenceToFileAndPath(os.path.join(karamba_plugin_path, "Karamba.gha"))
clr.AddReferenceToFileAndPath(os.path.join(karamba_plugin_path, "KarambaCommon.dll"))

import Karamba.Models.Model as Model
import Karamba.GHopper.Models.GH_Model as GH_Model
import Karamba.Elements.ModelTruss as Truss
import Karamba.Elements.ModelBeam as Beam
import Karamba.Materials.FemMaterial_Isotrop as FemMaterial
import feb # Karamba's C++ library (undocumented in the API)
import System.GC as GC

##############################

def randomize(iterable):
    sequence = list(iterable)
    random.shuffle(sequence)
    return sequence
    
def elapsed_time(start_time):
    return time.time() - start_time

##############################

def load_model(model):
    element_from_id = OrderedDict((e.ind, tuple(e.node_inds)) for e in model.elems)
    node_points = OrderedDict((n.ind, n.pos) for n in model.nodes)
    ground_nodes = {s.node_ind for s in model.supports}
    return element_from_id, node_points, ground_nodes

def get_id_from_element_from_model(model):
    return {tuple(e.node_inds) : e.ind for e in model.elems}

def nodes_from_elements(elements):
    return {n for e in elements for n in e}

def adjacent_from_edges(edges):
    """get a dict for graph connectivity {node_id : neighnor node ids}
    Parameters
    ----------
    edges : a list of 2-int-tuples
        [(node_1, node_2), ...]
    Returns
    -------
    undirected_edges : dict
        {node_id : neighnor node ids}
    """
    undirected_edges = defaultdict(set)
    for v1, v2 in edges:
        undirected_edges[v1].add(v2)
        undirected_edges[v2].add(v1)
    return undirected_edges

def get_distance(p1, p2):
    return math.sqrt(sum([(p1[i]-p2[i])**2 for i in range(3)]))

def incoming_from_edges(edges, element_from_id=None):
    incoming_vertices = defaultdict(set)
    for v1, v2 in edges:
        if element_from_id is None:
            incoming_vertices[v2].add(v1)
        else:
            incoming_vertices[element_from_id[v2]].add(element_from_id[v1])
    return incoming_vertices

##################################################

Node = namedtuple('Node', ['edge', 'vertex', 'cost'])

def dijkstra(source_vertices, successor_fn, cost_fn=lambda v1, v2: 1):
    node_from_vertex = {}
    queue = []
    for vertex in source_vertices:
        cost = 0
        node_from_vertex[vertex] = Node(None, None, cost)
        heapq.heappush(queue, (cost, vertex))
    while queue:
        cost1, vertex1 = heapq.heappop(queue)
        if node_from_vertex[vertex1].cost < cost1:
            continue
        for vertex2 in successor_fn(vertex1):
            edge = (vertex1, vertex2)
            cost2 = cost1 + cost_fn(*edge)
            if (vertex2 not in node_from_vertex) or (cost2 < node_from_vertex[vertex2].cost):
                node_from_vertex[vertex2] = Node(edge, vertex1, cost2)
                heapq.heappush(queue, (cost2, vertex2))
    return node_from_vertex

def compute_distance_from_node(elements, node_points, ground_nodes):
    nodes = nodes_from_elements(elements)
    neighbors = adjacent_from_edges(elements)
    edge_costs = {edge: get_distance(node_points[edge[0]], node_points[edge[1]])
                  for edge in elements}
    #edge_costs = {edge: np.abs(node_points[edge[1]] - node_points[edge[0]])[2]
    #              for edge in elements}
    edge_costs.update({edge[::-1]: distance for edge, distance in edge_costs.items()})
    successor_fn = lambda v: neighbors[v]
    cost_fn = lambda v1, v2: edge_costs[v1, v2] # TODO: z on cost function
    return dijkstra(ground_nodes & nodes, successor_fn, cost_fn)

#####################################

def retrace_sequence(visited, current_state, horizon=INF):
    (element, max_displacement), prev_state = visited[current_state]
    if (prev_state is None) or (horizon == 0):
        # tracing reaches the end
        return []
    previous_elements = retrace_sequence(visited, prev_state, horizon=horizon-1)
    return previous_elements + [(element, max_displacement)]

def compute_path_cost(plan):
    return sum([abs(max_disp) for _, max_disp in plan])

#####################################

def reverse_element(element):
    # element: int tuple
    return element[::-1]

def get_directions(element):
    return {element, reverse_element(element)}

def compute_candidate_elements(all_elements, ground_nodes, built_elements):
    nodes = compute_printed_nodes(ground_nodes, built_elements)
    for element in set(all_elements) - built_elements:
        for directed in get_directions(element):
            node1, node2 = directed
            if node1 in nodes:
                yield directed

def compute_printed_nodes(ground_nodes, printed):
    return nodes_from_elements(printed) | set(ground_nodes)

def nodes_from_elements(elements):
    return {n for e in elements for n in e}

def is_reversed(all_elements, element):
    assert (element in all_elements) != (reverse_element(element) in all_elements)
    return element not in all_elements

def get_undirected(all_elements, directed):
    is_reverse = is_reversed(all_elements, directed)
    assert (directed in all_elements) != is_reverse
    return reverse_element(directed) if is_reverse else directed

def get_id_from_element(element_from_id):
    return {e: i for i, e in element_from_id.items()}

def get_extructed_ids(element_from_id, directed_elements):
    id_from_element = get_id_from_element(element_from_id)
    extruded_ids = []
    for directed in directed_elements:
        element = get_undirected(id_from_element, directed)
        extruded_ids.append(id_from_element[element])
    return sorted(extruded_ids)

def get_other_node(node1, element):
    assert node1 in element
    return element[node1 == element[0]]

#######################################

def test_stiffness(model_in, element_from_id, elements, lc=0, verbose=False):
    if not elements:
        return True
    
#    # clone the model and its list of elements to avoid side effects
#    model = model_in.Clone()
#    # clone its elements to avoid side effects
#    model.cloneElements()
#    # clone the feb-model to avoid side effects
#    model.deepCloneFEModel()
    
    existing_e_ids = get_extructed_ids(element_from_id, elements)
    for e in model.elems:
        e.set_is_active(model, e.ind in existing_e_ids)
    
    deform = feb.Deform(model.febmodel)
    response = feb.Response(deform)
    
    try:
        # calculate the displacements
        response.updateNodalDisplacements();
        # calculate the member forces
        response.updateMemberForces();
    except:
        raise ValueError("The stiffness matrix of the system is singular.")
        
    # if something changed inform the feb-model about it (otherwise it won't recalculate)
    model.febmodel.touch()
    
#    energy_visitor = feb.EnergyVisitor(model.febmodel, model.febmodel.state(0), lc);
#    energy_visitor.visit(model.febmodel);
#    elastic_E = energy_visitor.elasticEnergy()
#    if verbose:
#        print '---'
##        print existing_e_ids
#        print 'nKinematicMode: ', deform.nKinematicModes()
#        print 'maxDisp (m): ', response.maxDisplacement()
#        print 'compliance: ', deform.compliance(lc)
#        print 'elastic E (kN-m): ', elastic_E
#        print '---'

#    axial_E = energy_visitor.axialEnergy()
#    bending_E = energy_visitor.bendingEnergy()
#    print 'axial + bending - elastic', abs(axial_E + bending_E - elastic_E)
    
    # this guards the objects from being freed prematurely
#    GC.KeepAlive(deform)
#    GC.KeepAlive(response)
    
    return response.maxDisplacement()

########################################

def check_connected(ground_nodes, built_elements):
    if not built_elements:
        return True
    node_neighbors = get_node_neighbors(built_elements)
    queue = deque(ground_nodes)
    visited_nodes = set(ground_nodes)
    visited_elements = set()
    while queue:
        node1 = queue.popleft()
        for element in node_neighbors[node1]:
            visited_elements.add(element)
            node2 = get_other_node(node1, element)
            if node2 not in visited_nodes:
                queue.append(node2)
                visited_nodes.add(node2)
    return built_elements <= visited_elements

def get_node_neighbors(elements):
    node_neighbors = defaultdict(set)
    for e in elements:
        n1, n2 = e
        node_neighbors[n1].add(e)
        node_neighbors[n2].add(e)
    return node_neighbors

def compute_z_distance(node_points, element):
    # Distance to a ground plane
    # Opposing gravitational force
    return get_midpoint(node_points, element)[2]

def get_midpoint(node_points, element):
    n1, n2 = element
    p1 = node_points[n1]
    p2 = node_points[n2]
    return [(p1[i]+p2[i])/2 for i in range(3)]

#######################################

SEARCH_HEURISTIC = {
    'random', 'z_distance', 'graph_distance',
}

def get_search_heuristic_fn(element_from_id, node_points, ground_nodes, heuristic='random', forward=True):
    sign = +1 if forward else -1
    all_elements = frozenset(element_from_id.values())
    distance_from_node = compute_distance_from_node(all_elements, node_points, ground_nodes)
    def h_fn(built_elements, directed):
        element = get_undirected(all_elements, directed)
        # return bias
        # lower bias will be dequed first
        if heuristic == 'random':
            return random.random()
        elif heuristic == 'z_distance':
            return sign*compute_z_distance(node_points, element)
        elif heuristic == 'graph_distance':
            node_costs = [distance_from_node[node].cost for node in element]
            return sign*(node_costs[0]+node_costs[1])/2
        raise NotImplementedError('Search heuristic ({}) not implemented, the only available ones are: {}'.format(
            heuristic, SEARCH_HEURISTIC))
    return h_fn

def add_successors(queue, all_elements, grounded_nodes, heuristic_fn, built_elements, incoming_from_element=None, verbose=False):
    remaining = all_elements - built_elements
    num_remaining = len(remaining) - 1
    assert 0 <= num_remaining, num_remaining
    candidate_elements = list(compute_candidate_elements(all_elements, grounded_nodes, built_elements))
    # TODO make incoming_from_element safe
#    if verbose: print('add successors: candidate elements: {}'.format(candidate_elements))
    for directed in randomize(candidate_elements):
        element = get_undirected(all_elements, directed)
        if not (incoming_from_element[element] <= built_elements):
            continue
        # compute bias value
        bias = heuristic_fn(built_elements, directed)
        # `num_remaining` can be seen as the search node's "inverse" depth in the search tree
        # heapq will prioritize poping ones with lower num_remaining value first, 
        # which means the node is on a deeper lever of the search tree
        priority = (num_remaining, bias, random.random())
        heapq.heappush(queue, (priority, built_elements, directed))

SearchState = namedtuple('SearchState', ['action', 'state'])

MAX_STATES_STORED = 20
RECORD_BT = True
RECORD_CONSTRAINT_VIOLATION = True

def progression_sequencing(model, heuristic='z_distance', stiffness=True, trans_tol=DEFAULT_TRANS_TOL, 
        partial_orders=None, verbose=False, backtrack_limit=INF, max_time=INF, optimize=True):
    start_time = time.time()
    
    ## * prepare data
    element_from_id, node_points, ground_nodes = load_model(model)
    id_from_element = get_id_from_element_from_model(model)
    all_elements = frozenset(element_from_id.values())
    assert len(ground_nodes) > 0
    heuristic_fn = get_search_heuristic_fn(element_from_id, node_points, ground_nodes, heuristic=heuristic)
    
    ## * initial state
    initial_built = frozenset()
    queue = []
    # visited keeps track of search history, maps () : SearchState
    visited = {initial_built : SearchState((None, None), None)}
    # TODO: check connectivity and stiffness for initial state
    # add candidate states from the initial state to the queue
    incoming_from_element = incoming_from_edges(partial_orders, element_from_id)
    add_successors(queue, all_elements, ground_nodes, heuristic_fn, initial_built, 
        incoming_from_element=incoming_from_element, verbose=verbose)
    
    plan = None
    min_remaining = len(element_from_id)
    num_evaluated = max_backtrack = stiffness_failures = 0
    
    #############################################
    bt_data = []  # backtrack history
    cons_data = []  # constraint violation history

    def snapshot_state(data_list=None, reason='', queue_log_cnt=0):
        """a lot of global parameters are used
        """
        cur_data = {}
        cur_data['num_evaluated'] = num_evaluated # iter
        cur_data['reason'] = reason
        cur_data['elapsed_time'] = elapsed_time(start_time)
        cur_data['min_remaining'] = min_remaining
        cur_data['max_backtrack'] = max_backtrack

        cur_data['stiffness_failures'] = stiffness_failures
        
        cur_data['backtrack'] = backtrack
        cur_data['total_q_len'] = len(queue)
        cur_data['chosen_element'] = element
        planned_elements = retrace_sequence(visited, built_elements)
        cur_data['planned_elements'] = planned_elements
        cur_data['queue'] = []
        if data_list is not None and len(data_list) < MAX_STATES_STORED:
            data_list.append(cur_data)
        return cur_data
    # end snapshot
    #############################################
    solutions = [] 
    while queue and (time.time() - start_time < max_time) :
        bias, built_elements, directed = heapq.heappop(queue)
        num_remaining = len(all_elements) - len(built_elements)
        element = get_undirected(all_elements, directed)
        
        num_evaluated += 1
        if num_remaining < min_remaining:
            min_remaining = num_remaining
        
        backtrack = num_remaining - min_remaining
        if backtrack > max_backtrack:
            max_backtrack = backtrack
            if RECORD_BT:
                print('max backtrack increased to {}'.format(max_backtrack))
                snapshot_state(bt_data, reason='Backtrack')
        if backtrack_limit < backtrack:
            print("Backtrack limit reached!")
            break # continue
        
        if verbose:
            print('Iteration: {} | Min Remain: {} | Built: {}/{} | Element: {} | Backtrack: {} | Stiffness failures: {} | Time: {:.3f}'.format(
                num_evaluated, min_remaining, len(built_elements), len(all_elements), element, backtrack, stiffness_failures, elapsed_time(start_time)))
        
        next_built_elements = built_elements | {element}
        
        # * check constraint
        if (next_built_elements in visited):
            if verbose: print('State visited before: {}'.format(next_built_elements))
            continue
        
        assert check_connected(ground_nodes, next_built_elements)
        max_displacement = INF
        if stiffness:
            max_displacement = test_stiffness(model, element_from_id, next_built_elements, verbose=verbose)
            if max_displacement > trans_tol:
                if verbose:
                    print("Stiffness constraint violated: max disp {:.4f} | tol {:.4f} [m]".format(max_displacement, trans_tol))
                # the stiffness constraint is violated
                stiffness_failures += 1
                if RECORD_CONSTRAINT_VIOLATION:
                    snapshot_state(cons_data, reason='stiffness_violation')
                continue
        
        # * record history
        # SearchState(action, current state)
        visited[next_built_elements] = SearchState((element, max_displacement), built_elements)
        if all_elements <= next_built_elements:
            min_remaining = 0
            plan = retrace_sequence(visited, next_built_elements)
            solutions.append(plan)
            if not optimize:
                break
            else:
                continue

        # * continue to the next search level, add candidates to queue
        add_successors(queue, all_elements, ground_nodes, heuristic_fn, next_built_elements, 
            incoming_from_element=incoming_from_element, verbose=verbose)
    
    data = {
        'planning_time' : elapsed_time(start_time),
        'num_evaluated': num_evaluated,
        'min_remaining': min_remaining,
        'max_backtrack': max_backtrack,
        'stiffness_failures': stiffness_failures,
        #
    }
    snapshot_data = {
        'backtrack_history': bt_data,
        'constraint_violation_history': cons_data,
    }
    solutions = sorted(solutions, key=lambda path: compute_path_cost(path))
    return solutions, data, snapshot_data


############################
# main
model = Model_in
# clone the model and its list of elements to avoid side effects
model = model.Clone()

_partial_orders = tree_to_list(partial_orders, lambda x:x) if partial_orders.DataCount > 0 else []
print(_partial_orders)

nodal_translation_tol = nodal_translation_tol or 5e-2 # meter
solutions, data, snapshot_data = progression_sequencing(model, heuristic=heuristic, trans_tol=nodal_translation_tol, 
    stiffness=check_stiffness, verbose=verbose, max_time=max_time, partial_orders=_partial_orders, optimize=optimize)

print(data)
print('Saved snapshot states (maximal stored state {}): BT state #{} | constraint violation state #{}'.format(
    MAX_STATES_STORED, len(snapshot_data['backtrack_history']), len(snapshot_data['constraint_violation_history'])))

id_from_element = get_id_from_element_from_model(model)
all_elements = frozenset(id_from_element.keys())
if len(solutions) != 0:
    print('{} plan found!'.format(len(solutions)))
    plan = solutions[0]
    plan = [[id_from_element[e], max_disp] for (e, max_disp) in plan]
    plan = list_to_tree(plan, source=[])
else:
    print('No plan found with time limit {}, displacement tol {:.4f} [m], please check out the snapshot states!'.format(
        max_time, nodal_translation_tol))

BT_state = DataTree[System.Int32]()
cons_violated_state = DataTree[System.Int32]()

bt_data = snapshot_data['backtrack_history']
cons_data = snapshot_data['constraint_violation_history']
for i, bdata in enumerate(bt_data):
    for e in bdata['planned_elements']:
        BT_state.Add(id_from_element[e[0]], GHPath(i))
#        BT_state.Add(e[1], GHPath(i, 1))
for i, cdata in enumerate(cons_data):
    for e in cdata['planned_elements']:
        cons_violated_state.Add(id_from_element[e[0]], GHPath(i))
#        cons_violated_state.Add(e[1], GHPath(i,1))
    cons_violated_state.Add(id_from_element[cdata['chosen_element']], GHPath(i))
