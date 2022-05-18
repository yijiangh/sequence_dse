import os
import logging
from collections import namedtuple, defaultdict, deque
import numpy as np
import pybullet_planning as pp
from pyconmech.frame_analysis import Node, Material

CLOSE_PT_TOL = 1e-5

#########################################

def get_logger(name):
    logger = logging.getLogger(name)

    try:
        from colorlog import ColoredFormatter
        formatter = ColoredFormatter("%(log_color)s%(levelname)-8s%(reset)s %(white)s%(message)s",
                                     datefmt=None,
                                     reset=True,
                                     log_colors={'DEBUG': 'cyan', 'INFO': 'green',
                                                 'WARNING': 'yellow',
                                                 'ERROR': 'red', 'CRITICAL': 'red',
                                                 }
                                     )
    except ImportError:
        formatter = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s | %(message)s')

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

    return logger

LOGGER = get_logger('seq_dse')

##########################################

# def name_from_id(id):
#     return 'E{}'.format(id)

# def id_from_name(name):
#     return int(name.split('E')[1])

####################################

def find_nodes(new_pts, cm_nodes, tol=CLOSE_PT_TOL):
    node_inds = list()
    for pt in new_pts:
        for n in cm_nodes:
            if pp.get_distance(pt, n.point) < tol:
                node_inds.append(n.node_ind)
                break
        else:
            # unseen node
            n_id = len(cm_nodes)
            cm_nodes.append(Node(pt, n_id, False)) # grounded flag
            node_inds.append(n_id)
            # print('Not found: #{} | {}'.format(n_id, pt))
    return cm_nodes, node_inds

def convex_combination(x1, x2, w=0.5):
    #assert 0 <= w <= 1
    return (1 - w) * x1 + (w * x2)

def closest_point_segments(line_point_1_1, line_point_1_2, line_point_2_1, line_point_2_2):
    from gurobipy import Model, GRB, GurobiError
    np_line1 = [np.array(line_point_1_1), np.array(line_point_1_2)]
    np_line2 = [np.array(line_point_2_1), np.array(line_point_2_2)]

    model = Model(name='qp_closest_points_between_segments')
    model.setParam(GRB.Param.OutputFlag, False)
    t1 = np.full(np_line1[0].shape, model.addVar(lb=0.0, ub=1.0, name="t1"))
    t2 = np.full(np_line2[0].shape, model.addVar(lb=0.0, ub=1.0, name="t2"))
    # difference = convex_combination(*np_line1, t1) - convex_combination(*np_line2, t2)
    difference = ((1 - t1) * np_line1[0] + t1 * np_line1[1]) - \
                 ((1 - t2) * np_line2[0] + t2 * np_line2[1])
    model.setObjective(sum(difference*difference), sense=GRB.MINIMIZE)
    try:
        model.optimize()
    except GurobiError as e:
        raise e
    return [convex_combination(*np_line1, t1[0].x), convex_combination(*np_line2, t2[0].x)]

def get_connector_from_elements(connectors, elements):
    # elements: list of indices
    connector_from_elements = defaultdict(set)
    for e in elements:
        for c in connectors:
            if e in c:
                connector_from_elements[e].add(c)
    return connector_from_elements

def get_element_neighbors(connectors, elements):
    # elements: list of indices
    connector_from_elements = get_connector_from_elements(connectors, elements)
    element_neighbors = defaultdict(set)
    for e in elements:
        for c in connector_from_elements[e]:
            element_neighbors[e].update(c)
        element_neighbors[e].remove(e)
    return element_neighbors

def check_connected(connectors, grounded_elements, printed_elements):
    """check if a given partial structure is connected to the ground
    Parameters
    ----------
    connectors : list of 2-int tuples
        each entry are the indices into the element set,
    grounded_elements : set
        grounded element ids
    printed_elements : set
        printed element ids
    Returns
    -------
    bool
        True if connected to ground
    """
    # TODO: for stability might need to check 2-connected
    if not printed_elements:
        return True
    printed_grounded_elements = set(grounded_elements) & set(printed_elements)
    if not printed_grounded_elements:
        return False
    element_neighbors = get_element_neighbors(connectors, printed_elements)
    queue = deque(printed_grounded_elements)
    visited_elements = set()
    while queue:
        n_element = queue.popleft()
        visited_elements.add(n_element)
        for element in element_neighbors[n_element]:
            if element in printed_elements and element not in visited_elements:
                visited_elements.add(element)
                queue.append(element)
    return printed_elements <= visited_elements

def find_grounded_elements(elements, grounded_nodes):
    grounded_elements = []
    for ename, element in elements.items():
        for pt in element.fem_data.axis_points:
            for grounded_node in grounded_nodes:
                if pp.get_distance(pt, grounded_node) < CLOSE_PT_TOL:
                    grounded_elements.append(ename)
    return grounded_elements

def find_connected_nodes(existing_elements, new_element, grounded_nodes):
    connected_nodes = []
    for element in existing_elements:
        for pt in element.fem_data.axis_points:
            for new_node in new_element.fem_data.axis_points:
                if pp.get_distance(pt, new_node) < CLOSE_PT_TOL:
                    for cp in connected_nodes:
                        if pp.get_distance(cp, new_node) < CLOSE_PT_TOL:
                            # duplicated
                            break
                    else:
                        connected_nodes.append(new_node)
    for new_node in new_element.fem_data.axis_points:
        for gn in grounded_nodes:
            if pp.get_distance(gn, new_node) < CLOSE_PT_TOL:
                for cp in connected_nodes:
                    if pp.get_distance(cp, new_node) < CLOSE_PT_TOL:
                        # duplicated
                        break
                else:
                    connected_nodes.append(new_node)
    return connected_nodes

######################################################

def get_midpoint(elements, e_id):
    return np.average([elements[e_id].fem_data.axis_points[i] for i in range(2)], axis=0)

def compute_z_distance(element_from_index, element):
    # Distance to a ground plane
    # Opposing gravitational force
    return get_midpoint(element_from_index, element)[2]

######################################################

def mkdir(path):
    if not os.path.exists(path):
        os.mkdir(path)
