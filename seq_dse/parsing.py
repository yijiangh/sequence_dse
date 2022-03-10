import os
import json
import numpy as np
from compas.geometry import is_colinear
from collections import namedtuple, defaultdict
from itertools import combinations

import pybullet_planning as pp
from pyconmech.frame_analysis import GravityLoad, Node, Element, Support, Material, CrossSec, Material, Joint, Model, LoadCase

from seq_dse.utils import CLOSE_PT_TOL
from seq_dse.data_structures import SeqElement, PlanningData, FEMData

HERE = os.path.dirname(__file__)
DATA_DIR = os.path.join(HERE, '..', 'data')

########################

def conmech_model_from_problem_name(problem_name):
    with open(os.path.join(DATA_DIR, problem_name, f'{problem_name}.json'), 'r') as f:
        data = json.load(f)
    return conmech_model_from_json_data(data)

def conmech_model_from_json_data(karamba_model_data):
    model = Model.from_data(karamba_model_data)
    loadcase = LoadCase.from_data(karamba_model_data['loadcases']['0'])
    return model

def get_fem_element_from_bar_id_from_model(model):
    return {e.elem_ind : [e.elem_ind] for e in model.elements}

def graph_from_conmech_model(model):
    cm_nodes = model.nodes
    elements = {e.elem_ind : SeqElement(name=e.elem_ind, fem_data=FEMData(axis_points=[
        cm_nodes[e.end_node_inds[i]].point for i in range(2)]), planning_data=None) \
        for e in model.elements}
    connectors_from_point_id = defaultdict(set)
    for e in model.elements:
        for v in e.end_node_inds:
            connectors_from_point_id[v].add(e.elem_ind)
    connectors = set()
    for v, v_connected_elements in connectors_from_point_id.items():
        for e1, e2 in combinations(list(v_connected_elements), 2):
            connectors.add(frozenset((e1, e2)))
    grounded_nodes = []
    for cm_node in cm_nodes:
        if cm_node.is_grounded:
            grounded_nodes.append(cm_node.point)
    return elements, list(connectors), grounded_nodes
