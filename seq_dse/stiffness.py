import heapq
import time
from itertools import combinations
from termcolor import cprint
from collections import namedtuple, defaultdict
import numpy as np
import random

import pybullet_planning as pp
from pyconmech.frame_analysis import StiffnessChecker
from pyconmech.frame_analysis import GravityLoad, Node, Element, Support, Material, CrossSec, Material, Joint, Model, LoadCase

from seq_dse.utils import check_connected, compute_z_distance, find_grounded_elements
from seq_dse.parsing import get_fem_element_from_bar_id_from_model

TRANS_TOL = 0.001 # 0.0005
ROT_TOL = np.inf # 5 * np.pi / 180

####################################

# TODO: get_max_nodal_deformation
Deformation = namedtuple('Deformation', ['success', 'displacements', 'fixities', 'reactions'])
Displacement = namedtuple('Displacement', ['dx', 'dy', 'dz', 'theta_x', 'theta_y', 'theta_z'])
Reaction = namedtuple('Reaction', ['fx', 'fy', 'fz', 'mx', 'my', 'mz'])

###################################################

def draw_conmech_model(model, fem_element_from_bar_id, existing_elements=None, exagg=10.0, deformations=None):
    # elements = elements or list(fem_element_from_bar_id.keys())
    existing_elements = list(range(len(model.elements))) if existing_elements is None else existing_elements
    cm_nodes = model.nodes
    cm_elements = model.elements
    cm_supports = model.supports
    h = []
    with pp.LockRenderer():
        existing_nodes = set()
        for e_id in existing_elements:
            for fem_e in fem_element_from_bar_id[e_id]:
                existing_nodes.update(set(cm_elements[fem_e].end_node_inds))
        # for n in cm_nodes:
        #     draw_point(n.point, size=0.002, color=GREEN if n.is_grounded else BLACK)
        for e_id in existing_elements:
            if e_id in fem_element_from_bar_id and len(fem_element_from_bar_id[e_id]) > 0:
                elem_ids = fem_element_from_bar_id[e_id]
            else:
                elem_ids = [e_id]
                # print('E#{}: {}'.format(e_id, elem_ids))
            for elem_id in elem_ids:
                e = cm_elements[elem_id]
                node_inds = e.end_node_inds
                if 'bar' in e.elem_tag or '' in e.elem_tag:
                    # (random.random(), random.random(), random.random(), 1)
                    th = pp.add_line(cm_nodes[node_inds[0]].point, cm_nodes[node_inds[1]].point,
                        color=pp.BLUE, width=3)
                    h.append(pp.add_text(e.elem_tag, position=(np.array(cm_nodes[node_inds[0]].point)+\
                        np.array(cm_nodes[node_inds[1]].point))/2))
                elif 'connector' in e.elem_tag:
                    th = pp.add_line(cm_nodes[node_inds[0]].point, cm_nodes[node_inds[1]].point, color=pp.BLUE, width=3)
                elif 'scaffold' in e.elem_tag:
                    node_inds = e.end_node_inds
                    th = pp.add_line(cm_nodes[node_inds[0]].point, cm_nodes[node_inds[1]].point, color=pp.TAN, width=3)
                    h.append(pp.add_text(e.elem_tag, position=(np.array(cm_nodes[node_inds[0]].point)+\
                        np.array(cm_nodes[node_inds[1]].point))/2))
                else:
                    raise ValueError("undefined elem tag {}".format(e.elem_tag))
                h.append(th)
                if deformations:
                    n1_d = deformations.displacements[node_inds[0]]
                    n2_d = deformations.displacements[node_inds[1]]
                    deform_h = pp.add_line(cm_nodes[node_inds[0]].point + exagg * np.array([n1_d.dx, n1_d.dy, n1_d.dz]),
                                     cm_nodes[node_inds[1]].point + exagg * np.array([n2_d.dx, n2_d.dy, n2_d.dz]),
                                     color=pp.BLACK, width=1)
                    h.append(deform_h)
        for _, s in cm_supports.items():
            if s.node_ind in existing_nodes:
                ph = pp.draw_point(cm_nodes[s.node_ind].point, size=0.2, color=pp.RED)
                h.extend(ph)
    return h

#########################################

# def create_stiffness_checker_from_graph(elements, connectors, grounded_nodes, verbose=False, trans_tol=TRANS_TOL, **kwargs):
#     # * convert structure to a conmech model
#     model, fem_element_from_bar_id = conmech_model_from_bar_structure(elements, connectors, grounded_nodes, **kwargs)
#     with pp.HideOutput():
#         checker = StiffnessChecker(model, verbose=verbose, checker_engine="numpy")
#         checker.set_loads(LoadCase(gravity_load=GravityLoad([0,0,-1.0])))
#     checker.set_nodal_displacement_tol(trans_tol=trans_tol, rot_tol=ROT_TOL)
#     return checker, fem_element_from_bar_id

def create_stiffness_checker_from_model(model, verbose=False, trans_tol=TRANS_TOL):
    with pp.HideOutput():
        checker = StiffnessChecker(model, verbose=verbose, checker_engine="numpy")
        checker.set_loads(LoadCase(gravity_load=GravityLoad([0,0,-1.0])))
    checker.set_nodal_displacement_tol(trans_tol=trans_tol, rot_tol=ROT_TOL)
    return checker, get_fem_element_from_bar_id_from_model(model)

##############################################

def evaluate_stiffness(model, existing_elements, checker=None, fem_element_from_bar_id=None, trans_tol=TRANS_TOL, verbose=True):
    # ! existing_elements uses int-based indexing
    # TODO: check each connected component individually
    # if not elements:
    #     return Deformation(True, {}, {}, {})

    if checker is None or fem_element_from_bar_id is None:
        checker, fem_element_from_bar_id = create_stiffness_checker_from_model(model, trans_tol=trans_tol, verbose=verbose)

    exist_element_ids = set()
    for elem in existing_elements:
        if elem in fem_element_from_bar_id:
            # ! bar
            exist_element_ids.update(fem_element_from_bar_id[elem])
        else:
            # ! in the case of scaffolding
            exist_element_ids.add(elem)

    is_stiff = checker.solve(exist_element_ids=list(exist_element_ids), if_cond_num=True)
    if not checker.has_stored_result():
        raise ValueError('checker no result')
        return Deformation(False, {}, {}, {})
    #print("has stored results: {0}".format(checker.has_stored_result()))

    success, nodal_displacement, fixities_reaction, element_reaction = checker.get_solved_results()
    assert is_stiff == success, "partial structure not stiff!"
    assert success or checker.get_compliance() > 0.0, 'partial structure: success {} | compliance {}'.format(success, checker.get_compliance())
    displacements = {i: Displacement(*d) for i, d in nodal_displacement.items()}
    fixities = {i: Reaction(*d) for i, d in fixities_reaction.items()}
    reactions = {i: (Reaction(*d[0]), Reaction(*d[1])) for i, d in element_reaction.items()}

    #print("nodal displacement (m/rad):\n{0}".format(nodal_displacement)) # nodes x 7
    #print("fixities reaction (kN, kN-m):\n{0}".format(fixities_reaction)) # ground x 7
    #print("element reaction (kN, kN-m):\n{0}".format(element_reaction)) # elements x 13
    trans_tol, rot_tol = checker.get_nodal_deformation_tol()
    max_trans, max_rot, max_trans_vid, max_rot_vid = checker.get_max_nodal_deformation()
    # The inverse of stiffness is flexibility or compliance

    translation = np.max(np.linalg.norm([d[:3] for d in displacements.values()], ord=2, axis=1))
    rotation = np.max(np.linalg.norm([d[3:] for d in displacements.values()], ord=2, axis=1))
    is_stiff &= (translation <= trans_tol) and (rotation <= rot_tol)

    if verbose:
        print('Stiff: {} | Compliance: {:.5f}'.format(is_stiff, checker.get_compliance()))
        print('Max translation deformation: {0:.5f} / {1:.5} = {2:.5}, at node #{3}'.format(
            max_trans, trans_tol, max_trans / trans_tol, max_trans_vid))
        # print('Max rotation deformation: {0:.5f} / {1:.5} = {2:.5}, at node #{3}'.format(
        #     max_rot, rot_tol, max_rot / rot_tol, max_rot_vid))

    return Deformation(is_stiff, displacements, fixities, reactions)

def test_stiffness(model, existing_elements, **kwargs):
    # ! existing_elements uses int-based indexing
    return evaluate_stiffness(model, existing_elements, **kwargs).success

##################################################

def plan_stiffness(model, elements, connectors, grounded_nodes, existing_elements, checker=None, fem_element_from_bar_id=None, trans_tol=TRANS_TOL, stiffness=True, heuristic='z', max_time=pp.INF, max_backtrack=0, verbose=False):
    """use the progression (foward-search) algorithm to plan a stiff sequence
    """
    start_time = time.time()
    # TODO the bar index gives the algorithm hints, try other starting point
    # TODO chosen bars
    if stiffness and (checker is None or fem_element_from_bar_id is None):
        checker, fem_element_from_bar_id = create_stiffness_checker_from_model(model, trans_tol=trans_tol)
    grounded_elements = find_grounded_elements(elements, grounded_nodes)

    # all_elements = frozenset(element_from_index)
    remaining_elements = frozenset(existing_elements)
    min_remaining = len(remaining_elements)
    max_bt = stiffness_failures = 0
    queue = [(None, frozenset(), [])]
    while queue and (pp.elapsed_time(start_time) < max_time):
        # TODO pop position and try distance heuristic
        _, printed, sequence = heapq.heappop(queue)
        num_remaining = len(remaining_elements) - len(printed)
        backtrack = num_remaining - min_remaining
        max_bt = max(max_bt, backtrack)
        if max_backtrack < backtrack:
            break # continue

        # * check constraints
        if not check_connected(connectors, grounded_elements, printed):
            continue
        if stiffness and not test_stiffness(elements, connectors, grounded_nodes, printed, checker=checker, fem_element_from_bar_id=fem_element_from_bar_id, verbose=verbose):
            stiffness_failures += 1
            continue

        if printed == remaining_elements:
            # * Done!
            #from extrusion.visualization import draw_ordered
            # distance = compute_sequence_distance(node_points, sequence, start=initial_position, end=initial_position)
            cprint('Plan-stiffness success! Elements: {}, max BT: {}, stiffness failure: {}, Time: {:.3f}sec'.format(
                len(sequence), max_bt, stiffness_failures, pp.elapsed_time(start_time))) #Distance: {:.3f}m,
            #local_search(extrusion_path, element_from_id, node_points, ground_nodes, checker, sequence,
            #             initial_position=initial_position, stiffness=stiffness, max_time=INF)
            #draw_ordered(sequence, node_points)
            #wait_for_user()
            return sequence

        # * add successors
        for element in pp.randomize(remaining_elements - printed):
            new_printed = printed | {element}
            new_sequence = sequence + [element]
            num_remaining = len(remaining_elements) - len(new_printed)
            min_remaining = min(min_remaining, num_remaining)
            # Don't count edge length
            # distance = get_distance(position, node_points[node1]) if position is not None else None
            # distance = compute_sequence_distance(node_points, new_sequence)
            if heuristic == 'none':
                bias = None
            elif heuristic == 'random':
                bias = random.random()
            elif heuristic == 'z':
                bias = compute_z_distance(elements, element)
            # elif heuristic == 'distance':
            #     bias = distance
            else:
                raise ValueError(heuristic)
            #bias = heuristic_fn(printed, element, conf=None) # TODO: experiment with other biases
            priority = (num_remaining, bias, random.random())
            heapq.heappush(queue, (priority, new_printed, new_sequence))

    cprint('Failed to find stiffness plan under tol {}! Elements: {}, Min remaining {}, Time: {:.3f}sec'.format(
        trans_tol, len(remaining_elements), min_remaining, pp.elapsed_time(start_time)), 'red')
    return None

#########################################

def compute_plan_deformation(model, plan, verbose=False, trans_tol=TRANS_TOL, checker=None,
    fem_element_from_bar_id=None):
    if checker is None or fem_element_from_bar_id is None:
        checker, fem_element_from_bar_id = create_stiffness_checker_from_model(model, verbose=False, trans_tol=trans_tol)
    trans_tol, rot_tol = checker.get_nodal_deformation_tol()
    if plan is None:
        return trans_tol, rot_tol

    printed = []
    translations = []
    rotations = []
    for i, element in enumerate(plan):
        if verbose:
            cprint('{}) : bar#{}'.format(i, element), 'yellow')
        printed.append(element)
        deformation = evaluate_stiffness(model, printed,
                                         checker=checker, fem_element_from_bar_id=fem_element_from_bar_id, verbose=verbose)
        trans, rot, _, _ = checker.get_max_nodal_deformation()
        translations.append([bool(deformation.success), trans])
        rotations.append(rot)
    # TODO: could return full history
    return translations
