import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

import pybullet_planning as pp
# from utils import 
from seq_dse.stiffness import evaluate_stiffness, draw_conmech_model
from seq_dse.parsing import conmech_model_from_problem_name, graph_from_conmech_model, get_fem_element_from_bar_id_from_model

#########################################

def main():
    np.set_printoptions(precision=6)
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--viewer', action='store_true', help='Enables the viewer during planning (slow!)')
    parser.add_argument('-p', '--problem', default='roof', help='Problem name to be solved.')
    #
    parser.add_argument('--debug', action='store_true', help='Debug verbose mode')

    args = parser.parse_args()
    print('Arguments:', args)

    pp.connect(use_gui=args.viewer)
    camera_base_pt = (0, 0, 1)
    camera_pt = np.array(camera_base_pt) + np.array([0, 2, 1.5])
    pp.set_camera_pose(tuple(camera_pt), camera_base_pt)
    pp.draw_pose(pp.unit_pose())

    model = conmech_model_from_problem_name(args.problem)
    fem_element_from_bar_id = get_fem_element_from_bar_id_from_model(model)
    # conmech_model_from_bar_structure(elements, connectors, grounded_nodes, debug=True, save_model_path=None)
    elements, connectors, grounded_nodes = graph_from_conmech_model(model)

    if pp.has_gui() and args.debug:
        existing_elements = list(elements)
        deformations = evaluate_stiffness(model, existing_elements,
            verbose=True)
        # deformations = None
        h = draw_conmech_model(model, fem_element_from_bar_id, deformations=deformations, exagg=100)
        pp.wait_if_gui()
        with pp.LockRenderer():
            pp.remove_handles(h)

    # sequence = plan_stiffness(elements, connectors, grounded_nodes, existing_elements)
    # wait_if_gui('Ready to simulate construction sequence.')
    # if pp.has_gui():
    #     remove_all_debug()
    #     # checker, fem_element_from_bar_id = create_stiffness_checker(elements, connectors, grounded_nodes)
    #     for i in range(len(sequence)):
    #         pp.remove_all_debug()
    #         h = draw_conmech_model(model, fem_element_from_bar_id, existing_elements=sequence[:i])
    #         connected_nodes = find_connected_nodes([elements[name_from_id(te)] for te in sequence[:i]], elements[name_from_id(sequence[i])], grounded_nodes)
    #         for cn in connected_nodes:
    #             pp.draw_point(cn, color=pp.GREEN, size=0.2)
    #         print('Step {}:'.format(i))
    #         wait_if_gui()
            # evaluate_stiffness(model, sequence[:i], 
                # checker=checker, fem_element_from_bar_id=fem_element_from_bar_id, verbose=True)
            # wait_for_duration(0.0)

    # translations = compute_plan_deformation(elements, connectors, grounded_nodes, sequence)
    # plt.figure(1)
    # fig, axs = plt.subplots(1,2)
    # axs[0].plot([t[1] for t in translations])
    # axs[0].plot([TRANS_TOL for _ in translations])
    # axs[0].set_xlabel('Construction Iter')
    # axs[0].set_ylabel('Max nodal deformation (m)')
    # axs[0].set_title('Progression deformation history')

    # axs[1].plot(sequence)
    # axs[1].set_xlabel('Construction Iter')
    # axs[1].set_ylabel('Element index')
    # axs[1].set_title('sequence')
    # fig.show()

    pp.wait_for_user('Done.')
    pp.reset_simulation()
    pp.disconnect()

if __name__ == '__main__':
    main()