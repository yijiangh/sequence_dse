import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

import pybullet_planning as pp
from seq_dse.stiffness import evaluate_stiffness, draw_conmech_model, plan_stiffness, compute_plan_deformation, TRANS_TOL
from seq_dse.parsing import conmech_model_from_problem_name, graph_from_conmech_model, get_fem_element_from_bar_id_from_model
from seq_dse.utils import find_connected_nodes

#########################################

def main():
    np.set_printoptions(precision=6)
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--viewer', action='store_true', help='Enables the viewer during planning (slow!)')
    parser.add_argument('-p', '--problem', default='roof', help='Problem name to be solved.')
    parser.add_argument('--optimize', action='store_true', help='Find the optimal sequence.')
    #
    parser.add_argument('--debug', action='store_true', help='Debug verbose mode')

    args = parser.parse_args()
    print('Arguments:', args)

    pp.connect(use_gui=args.viewer)
    camera_base_pt = (0, 0, 1)
    camera_pt = np.array(camera_base_pt) + np.array([0, 20, 1.5])
    pp.set_camera_pose(tuple(camera_pt), camera_base_pt)
    pp.draw_pose(pp.unit_pose())

    model, loadcase = conmech_model_from_problem_name(args.problem)
    fem_element_from_bar_id = get_fem_element_from_bar_id_from_model(model)
    # conmech_model_from_bar_structure(elements, connectors, grounded_nodes, debug=True, save_model_path=None)
    elements, connectors, grounded_nodes = graph_from_conmech_model(model)
    starting_existing_elements = list(elements)

    if pp.has_gui() and args.debug:
        deformations = evaluate_stiffness(model, starting_existing_elements,
            verbose=True) #, loadcase=loadcase)
        print('Nodal displacements')
        for i, reaction in deformations.displacements.items():
            print(f'N{i} : {reaction.dx} | {reaction.dy} | {reaction.dz}')
        print('Element axial force')
        for i, reaction in deformations.reactions.items():
            print(f'E{i} : {reaction[0].fx} | {reaction[1].fx}')
        # deformations = None
        h = draw_conmech_model(model, fem_element_from_bar_id, deformations=deformations, exagg=100)
        pp.wait_if_gui()
        with pp.LockRenderer():
            pp.remove_handles(h)
        pp.disconnect()
        return

    sequence = plan_stiffness(model, elements, connectors, grounded_nodes, starting_existing_elements, verbose=True)
    if sequence:
        pp.wait_if_gui('Ready to simulate construction sequence.')
        if pp.has_gui():
            pp.remove_all_debug()
            # checker, fem_element_from_bar_id = create_stiffness_checker(elements, connectors, grounded_nodes)
            for i in range(len(sequence)):
                with pp.LockRenderer():
                    pp.remove_all_debug()
                    h = draw_conmech_model(model, fem_element_from_bar_id, existing_elements=sequence[:i])
                    # * draw connected nodes of the current element
                    connected_nodes = find_connected_nodes([elements[te] for te in sequence[:i]], elements[sequence[i]], grounded_nodes)
                    for cn in connected_nodes:
                        pp.draw_point(cn, color=pp.GREEN, size=0.2)
                print('Step {}:'.format(i))
                # evaluate_stiffness(model, sequence[:i], 
                #     checker=checker, fem_element_from_bar_id=fem_element_from_bar_id, verbose=True)
                pp.wait_if_gui()
                pp.wait_for_duration(0.05)

        translations = compute_plan_deformation(model, sequence)
        plt.figure(1)
        fig, axs = plt.subplots(1,2)
        axs[0].plot([t[1] for t in translations])
        axs[0].plot([TRANS_TOL for _ in translations])
        axs[0].set_xlabel('Construction Iter')
        axs[0].set_ylabel('Max nodal deformation (m)')
        axs[0].set_title('Progression deformation history')

        axs[1].plot(sequence)
        axs[1].set_xlabel('Construction Iter')
        axs[1].set_ylabel('Element index')
        axs[1].set_title('sequence')
        fig.show()

    pp.wait_for_user('Done.')
    pp.reset_simulation()
    pp.disconnect()

if __name__ == '__main__':
    main()