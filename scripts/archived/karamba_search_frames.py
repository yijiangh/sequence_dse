# The much of the code in the sequencing part is adapted from https://github.com/caelan/pb-construction
if not enable:
    raise ValueError('Search not enabled.')

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

def randomize(iterable):
    sequence = list(iterable)
    random.shuffle(sequence)
    return sequence
    
def elapsed_time(start_time):
    return time.time() - start_time

##############################

def load_model(model, grounded_supports):
    # load connectivity data from karamba model
    element_from_id = OrderedDict((e.ind, tuple(e.node_inds)) for e in model.elems)
    node_points = OrderedDict((n.ind, n.pos) for n in model.nodes)
#    ground_nodes = {s.node_ind for s in model.supports}
    ground_nodes = set()
    for k_supp in model.supports:
        for k_g in grounded_supports:
            if k_supp.position.DistanceTo(k_g.position) < 1e-6:
              ground_nodes.add(k_supp.node_ind)
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

def get_sequence_from_plan(plan):
    return [p[0] for p in plan]

#####################################

def reverse_element(element):
    # element: int tuple
    return element[::-1]

def get_directions(element):
    return {element, reverse_element(element)}

def compute_candidate_elements(all_elements, ground_nodes, built_elements, use_two_link_prune=use_two_link_prune):
    # connectivity constraint is enforced here
    # enforce two link constraint
    elements_from_node = get_node_neighbors(built_elements)
    nodes = compute_printed_nodes(ground_nodes, built_elements)
    for element in set(all_elements) - built_elements:
#        for directed in get_directions(element):
#        directed = element
        node1, node2 = element
        if (node1 not in nodes) and (node2 not in nodes):
            # element floating
            continue
        if not use_two_link_prune:
            yield element
        else:
            if not (node1 in nodes and node2 in nodes):
                if node1 in nodes:
                    # if existing node is not grounded and is degree-one before adding the new element
                    # aka a two-link cantilever
                    if node1 not in ground_nodes and len(elements_from_node[node1]) == 1:
                        continue
                else:
                    if node2 in ground_nodes and len(elements_from_node[node2]) == 1:
                        continue
            yield element

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

########################################

def check_connected(ground_nodes, seq_built_elements):
    if not seq_built_elements:
        return True
    built_elements = set(seq_built_elements)
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
    def h_fn(built_elements, element):
#        element = get_undirected(all_elements, directed)
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

def add_successors(queue, all_elements, grounded_nodes, heuristic_fn, ordered_built_elements, 
        incoming_from_element=None, verbose=False, solutions=[], violated_states=[]):
    built_elements = set(ordered_built_elements)
    remaining = all_elements - built_elements
    num_remaining = len(remaining) - 1
    assert 0 <= num_remaining, 'num_remaining {}'.format(num_remaining)
    candidate_elements = list(compute_candidate_elements(all_elements, grounded_nodes, built_elements))
    # TODO make incoming_from_element safe
    for directed in candidate_elements: #randomize(candidate_elements):
        element = get_undirected(all_elements, directed)
        # partial ordering
        if not (incoming_from_element[element] <= built_elements):
            continue
        # compute bias value
        if len(solutions) == 0:
            bias = heuristic_fn(built_elements, element)
        else:
            # deflation-like objective
            new_ordered_built_elements = ordered_built_elements + (element,)
            similarity_score = 0
            for _, prev_plan in solutions:
                prev_seq = get_sequence_from_plan(prev_plan)
                for q_i, q_e in enumerate(new_ordered_built_elements):
                    similarity_score += abs(q_i - prev_seq.index(q_e))
            similarity_score *= -1.0
            bias = similarity_score
        # `num_remaining` can be seen as the search node's "inverse" depth in the search tree
        # heapq will prioritize poping ones with lower num_remaining value first, 
        # which means the node is on a deeper lever of the search tree
        priority = (num_remaining, bias) #, random.random())
        heapq.heappush(queue, (0, priority, ordered_built_elements, element))

SearchState = namedtuple('SearchState', ['action', 'state'])

#MAX_STATES_STORED = 20
#RECORD_BT = True
#RECORD_CONSTRAINT_VIOLATION = False

def progression_sequencing(model, grounded_supports, heuristic='z_distance', stiffness=True, trans_tol=DEFAULT_TRANS_TOL, 
        partial_orders=None, verbose=False, backtrack_limit=INF, max_time=INF, optimize=True):
    start_time = time.time()
    
    ## * prepare data
    element_from_id, node_points, ground_nodes = load_model(model, grounded_supports)
    id_from_element = get_id_from_element_from_model(model)
    all_elements = frozenset(element_from_id.values())
    assert len(ground_nodes) > 0
    heuristic_fn = get_search_heuristic_fn(element_from_id, node_points, ground_nodes, heuristic=heuristic)
    
    ## * initial state
    initial_built = ()
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
    cached_analysis = {}
    
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

    best_solution_score = 1e8
    solutions = []
    while queue and (time.time() - start_time < max_time) :
        level, bias, built_elements, element = heapq.heappop(queue)
        num_remaining = len(all_elements) - len(built_elements)
#        element = get_undirected(all_elements, directed)
        
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
            print('Iteration: {} | Min Remain: {} | Built: {}/{} | Level: {} | PB: {} | Element: {} | Backtrack: {} | Stiffness failures: {} | Time: {:.3f}'.format(
                num_evaluated, min_remaining, len(built_elements), len(all_elements), 
                level, [id_from_element[e] for e in built_elements], id_from_element[element], backtrack, stiffness_failures, elapsed_time(start_time)))
        
#        next_built_elements = built_elements | {element}
        next_built_elements = built_elements + (element,)
        
        # * check constraint
        if len(solutions) > 0 and (next_built_elements in visited):
            if verbose: print('State visited before: {}'.format(next_built_elements))
            continue
        
#        assert check_connected(ground_nodes, next_built_elements)
        max_displacement = INF
        if stiffness:
            built_element_set = frozenset(next_built_elements)
            if built_element_set not in cached_analysis:
                existing_e_ids = get_extructed_ids(element_from_id, next_built_elements)
                max_displacement = test_stiffness_fn(model, existing_e_ids, verbose=verbose)
                cached_analysis[built_element_set] = max_displacement
            else:
                if verbose:
                    print('reuse analysis!')
                max_displacement = cached_analysis[built_element_set]
            if max_displacement > trans_tol:
                if verbose:
                    print("Stiffness constraint violated: max val {:.4f} | tol {:.4f}".format(
                    max_displacement, trans_tol))
                # the stiffness constraint is violated
                stiffness_failures += 1
                if RECORD_CONSTRAINT_VIOLATION:
                    snapshot_state(cons_data, reason='stiffness_violation')
                continue
        
            # check if max_displacement is bigger than the best score, we skip
            # TODO generalize to arbitrary path_cost_fn
            if optimize and (not store_worse_solution) and max_displacement > best_solution_score:
#                print('Prune worse-performing branch {:.4f} | best {:.4f}.'.format(max_displacement, best_solution_score))
                continue
        
        # * record history
        # SearchState(action, current state)
        visited[next_built_elements] = SearchState((element, max_displacement), built_elements)
        if all_elements <= set(next_built_elements):
            min_remaining = 0
            # a solution has been found!
            plan = retrace_sequence(visited, next_built_elements)
            seq = get_sequence_from_plan(plan)
            new_path_cost = compute_path_cost_fn(plan)
            if not optimize:
                heapq.heappush(solutions, (new_path_cost, plan))
                break
            else:
                # check duplicated plan
                duplicated_plan = False
                for _, prev_plan in solutions:
                    if get_sequence_from_plan(prev_plan) == seq:
                        print('> Plan found before: {}'.format([id_from_element[e] for e in seq]))
                        duplicated_plan = True
                        break
                # exit the search if a duplicated plan is found
                if duplicated_plan:
                    break
                if store_worse_solution or new_path_cost <= best_solution_score*(1+1e-6):
                    print('!!! score {:.6f} | #of sols: {}|  better plan found: {}'.format(new_path_cost, len(solutions), 
                        [id_from_element[e] for e in seq]))
                    best_solution_score = new_path_cost
                    # TODO can also not append if not a better solution
                    heapq.heappush(solutions, (new_path_cost, plan))
                else:
                    print('? score {:.6f} | #of sols: {}| not a better plan found: {}'.format(new_path_cost, len(solutions), 
                        [id_from_element[e] for e in seq]))
                # clear the queue and restart the search
                queue = []
                visited = {initial_built : SearchState((None, None), None)}
                add_successors(queue, all_elements, ground_nodes, heuristic_fn, initial_built, 
                    incoming_from_element=incoming_from_element, verbose=verbose, solutions=solutions)
                continue
                
        # * continue to the next search level, add candidates to queue
        add_successors(queue, all_elements, ground_nodes, heuristic_fn, next_built_elements, 
            incoming_from_element=incoming_from_element, verbose=verbose, solutions=solutions)
    
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
#    solutions = sorted(solutions, key=lambda path: compute_path_cost_fn(path))
    return solutions, data, snapshot_data


############################
# main

if enable:
    model = Model_in
    # clone the model and its list of elements to avoid side effects
    model = model.Clone()
    
    _partial_orders = tree_to_list(partial_orders, lambda x:x) if partial_orders.DataCount > 0 else []
#    print(_partial_orders)
    
#    nodal_translation_tol = nodal_translation_tol or 5e-2 # meter
    element_solutions, data, snapshot_data = progression_sequencing(model, grounded_supports, 
        heuristic=heuristic, trans_tol=nodal_translation_tol,
        stiffness=check_stiffness, verbose=verbose, max_time=max_time, partial_orders=_partial_orders, optimize=optimize)
    
    print(data)
    print('Saved snapshot states (maximal stored state {}): BT state #{} | constraint violation state #{}'.format(
        MAX_STATES_STORED, len(snapshot_data['backtrack_history']), len(snapshot_data['constraint_violation_history'])))
    
    id_from_element = get_id_from_element_from_model(model)
    all_elements = frozenset(id_from_element.keys())
    
    solutions = [[(id_from_element[e], d) for e, d in e_path] for score, e_path in element_solutions]
    solution_costs = [sdata[0] for sdata in element_solutions]

    if len(solutions) != 0:
        print('{} plan found!'.format(len(solutions))) 
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