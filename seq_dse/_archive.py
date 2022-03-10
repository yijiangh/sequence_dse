####################################
# E, G12, fy, density
# ! NOT USED NOW: fy is the material strength in the specified direction (local x direction)
# element_tags = None is the fall-back material entry

def wood_material(elem_tags=None):
    return Material(1050*1e4, 360*1e4, 1.3*1e4, 6, elem_tags=elem_tags, family='Wood', name='Wood', type_name='ISO')
# STEEL_MATERIAL = Material(1050*1e4, 360*1e4, 1.3*1e4, 78.5, elem_tags=None, family='Steel', name='Steel', type_name='ISO')

CONNECTOR_STENGTH_RATIO = 0.1
def connector_material(elem_tags=None):
    return Material(1050*1e4*CONNECTOR_STENGTH_RATIO, 360*1e4*CONNECTOR_STENGTH_RATIO, 1.3*1e4, 6,
        elem_tags=elem_tags, family='Connector', name='Glue', type_name='ISO')

ROTATIONAL_STIFFNESS = 100 # kn/rad

def A_solid_cir(radius):
    return np.pi * radius**2

def Jx_solid_cir(radius):
    # https://en.wikipedia.org/wiki/Polar_moment_of_inertia
    return np.pi * radius**4 / 2

def Iy_solid_cir(radius):
    # https://www.engineeringtoolbox.com/area-moment-inertia-d_1328.html
    return np.pi * radius**4 / 4

def solid_cir_crosssec(r, elem_tags=None):
    A = A_solid_cir(r)
    Jx = Jx_solid_cir(r)
    Iy = Iy_solid_cir(r)
    return CrossSec(A, Jx, Iy, Iy, elem_tags=elem_tags, family='Circular', name='solid_circle')

#########################################

def conmech_model_from_bar_structure(elements, connectors, grounded_nodes, debug=False, save_model_path=None):
    unit = "meter"
    # * Nodes
    cm_nodes = []
    cm_elements = []
    pts_from_element = defaultdict(list)
    fem_element_from_bar_id = defaultdict(set)
    for e_name, element in elements.items():
        e_id = id_from_name(e_name)
        # Node(self, point, node_ind, is_grounded):
        # ! is_grounded is used to mark grounded nodes for construction purpose
        # structural supports are specified using the Support entries
        cm_nodes.append(Node(element.fem_data.axis_points[0], e_id*2, False))
        cm_nodes.append(Node(element.fem_data.axis_points[1], e_id*2+1, False))
        pts_from_element[e_id].extend(element.fem_data.axis_points)
    # original_bar_cm_nodes = copy(cm_nodes)

    # * Supports
    supports = []
    cn_cnt = len(cm_nodes)
    cm_nodes, grounded_node_inds = find_nodes(grounded_nodes, cm_nodes)
    assert len(cm_nodes) == cn_cnt, 'supports not well defined!'
    for gn_id in grounded_node_inds:
        # Support(self, condition, node_ind)
        supports.append(Support([1 for _ in range(6)], gn_id))

    # * contact connectors
    contact_from_connectors = {}
    for c in connectors:
        e1 = elements[name_from_id(c[0])]
        e2 = elements[name_from_id(c[1])]
        connect_segment = closest_point_segments(*e1.fem_data.axis_points, *e2.fem_data.axis_points)
        # if all the bars are connected at their ends, there should be no gap between the lines
        if pp.get_distance(*connect_segment) > CLOSE_PT_TOL:
            contact_from_connectors[c] = connect_segment

    connector_tags = set()
    for c_id, connector_endpts in contact_from_connectors.items():
        if pp.get_distance(*connector_endpts) < CLOSE_PT_TOL:
            # no extra connector element needed
            continue
        cm_nodes, node_inds = find_nodes(connector_endpts, cm_nodes)
        for pt in connector_endpts:
            for e_id in c_id:
                # e_id != GROUND_INDEX and 
                if is_colinear(pt, *elements[name_from_id(e_id)].fem_data.axis_points, tol=CLOSE_PT_TOL):
                    pts_from_element[e_id].append(pt)
        assert len(node_inds) == 2
        assert node_inds[0] != node_inds[1], "#{}|{} - {}".format(c_id, connector_endpts, node_inds)
        # * add connector element
        # Element(self, end_node_inds, elem_ind, elem_tag='', bending_stiff=True)
        e_id = len(cm_elements)
        e_tag = 'connector'# .format(e_id)
        cm_elements.append(Element(tuple(node_inds), e_id, elem_tag=e_tag, bending_stiff=True))
        connector_tags.add(e_tag)
        for e in c_id:
            fem_element_from_bar_id[e].add(e_id)

    bar_tags = set()
    for e_id, seg_pts in pts_from_element.items():
        # assert len(seg_pts) == 2, 'e_id {} pts {}'.format(e_id, seg_pts)
        sorted_pt_pairs = sorted(combinations(list(seg_pts), 2), key=lambda pt_pair: pp.get_distance(*pt_pair))
        e_end_pts = elements[name_from_id(e_id)].fem_data.axis_points
        farthest_pts = sorted_pt_pairs[-1]
        if pp.get_angle(pp.get_difference(*farthest_pts), pp.get_difference(*e_end_pts)) > np.pi/2:
            farthest_pts = farthest_pts[::-1]
        sorted_seg_pts = sorted(list(seg_pts), key=lambda pt: pp.get_distance(pt, farthest_pts[0]))

        cm_nodes, node_inds = find_nodes(list(sorted_seg_pts), cm_nodes)
        assert len(sorted_seg_pts) == len(node_inds)
        elem_tag = 'bar' #.format(e_id)
        bar_tags.add(elem_tag)
        seg_id = 0
        for i in range(len(sorted_seg_pts)-1):
            if node_inds[i] != node_inds[i+1]:
                new_element = Element((node_inds[i], node_inds[i+1]), len(cm_elements), elem_tag=elem_tag, bending_stiff=True)
                cm_elements.append(new_element)
                fem_element_from_bar_id[e_id].add(new_element.elem_ind)
                seg_id += 1

    # TODO Add rotational stiffness later
    # ? Joint(self, c_conditions, elem_tags):
    node_c_conditions = [None, None, None] + [ROTATIONAL_STIFFNESS for _ in range(3)]
    # joint = Joint(node_c_conditions+node_c_conditions, list(bar_tags | connector_tags))
    joints = []

    # TODO different material property and cross secs on Element and Connectors
    r = elements[name_from_id(0)].fem_data.radius # in meter
    crosssecs = [solid_cir_crosssec(r, elem_tags=list(bar_tags | connector_tags))]
    materials = [wood_material(list(bar_tags)), connector_material(list(connector_tags))]

    model = Model(cm_nodes, cm_elements, supports, joints, materials, crosssecs, unit=unit)
    if save_model_path is not None:
        model_path = os.path.join(save_model_path, model.model_name.split(".json")[0] + '_conmech_model.json')
        model_data = model.to_data()
        model_data['fem_element_from_bar_id'] = {bar : list(fem_es) for bar, fem_es in fem_element_from_bar_id.items()}
        with open(model_path, 'w') as f:
            json.dump(model_data, f, indent=None)
        cprint('Conmech model saved to: {}'.format(model_path), 'green')

    return model, fem_element_from_bar_id

