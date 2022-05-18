import os
import json
import colorsys
import numpy as np
import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt
from .utils import mkdir

def sample_colors(num, lower=0.0, upper=1.0): # for now wrap around
    # lower=0.0, upper=0.75
    # return [colorsys.hsv_to_rgb(h, s=1, v=1) for h in reversed(np.linspace(lower, upper, num, endpoint=True))]
    viridis = cm.get_cmap('viridis', 12)
    return viridis(np.linspace(lower, upper, num, endpoint=True))

def save_seq_performance_plot(data_path: str, test_folder_name: str, design_id: int, seq_id: int, 
        seq, stress_history, disp_history):
    # csv_file_name = os.path.join(data_path, 'samples.csv')
    # orig_pd_data = pd.read_csv(csv_file_name, header=None)
    # orig_pd_data.columns = ['var1', 'var2', 'var3', 'obj', 'num_elem', 'weight', 'karamba_model_path', 'seq_data_path']
    # seq_datapath  = orig_pd_data["seq_data_path"][design_id]
    # with open(os.path.join(data_path, seq_datapath), 'r') as fp:
    #     seq_data = json.load(fp) 
    # seq = sorted(seq_data['solutions'], key=lambda x:x['seq_cost'])[seq_id]['seq']

    given_constraint_type = test_folder_name.split('_')[0].split('=')[0].split('-')[0]
    histories = {
        'displacement' : disp_history, 
        'stress' : stress_history,
    }

    pdf_path = os.path.join(data_path, 'raw_imgs')
    mkdir(pdf_path)
    for constraint_type, fea_history in histories.items():
        unit = '[kN/cm2]' if constraint_type == 'stress' else '[cm]'
        scale = 1.0 if constraint_type == 'stress' else 100.0

        plt.rc('font', size=14)
        # plt.clf()
        fig, ax = plt.subplots(figsize=(10, 2.5), layout='constrained')

        fea_history = scale * np.array(fea_history)
        sorted_ids = np.argsort(-fea_history)
        max_id = sorted_ids[0]
        max_fea = fea_history[max_id]

        # * performance line
        perf_color = '#E39774' if constraint_type == 'displacement' else '#326273'
        ax.plot(range(len(seq)), fea_history, c=perf_color, linewidth = 1)

        # * tolerance line
        if constraint_type == given_constraint_type:
            tol = float(test_folder_name.split('_')[0].split('=')[-1])
            ax.plot(range(len(seq)), [tol for step in seq], c='black', linestyle='--', linewidth = 0.5)
            # ax.plot(range(len(seq)), [max_fea for step in seq], c='black', linestyle='-', linewidth = 0.5)

        ax.set_xlabel('assembly steps')
        ax.set_ylabel('maximal {} {}'.format(constraint_type, unit))

        plt.annotate('{:.2f}'.format(max_fea), xy=(max_id, max_fea), xytext=(max_id, max_fea+2), xycoords='data',
                     arrowprops=dict(arrowstyle="->", lw=0.5),
                     )
        # extraticks = [max_fea]
        # plt.yticks(list(plt.yticks()[0]) + extraticks);

        fig.savefig(os.path.join(pdf_path, 'design{}_seq{}_{}.png'.format(design_id, seq_id, constraint_type)), dpi=200)
        fig.savefig(os.path.join(pdf_path, 'design{}_seq{}_{}.pdf'.format(design_id, seq_id, constraint_type)), dpi=200)

    return True