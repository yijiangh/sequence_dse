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
        seq, stress_history, disp_history, black_background=False):
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

    plt.rcParams['font.sans-serif'] = "Arial"
    plt.rcParams['font.family'] = "sans-serif"
    plt.rc('font', size=14)
    if black_background:
        plt.rcParams['axes.facecolor'] = 'black'
        COLOR = 'white'
        plt.rcParams['text.color'] = COLOR
        plt.rcParams['axes.labelcolor'] = COLOR
        plt.rcParams['xtick.color'] = COLOR
        plt.rcParams['ytick.color'] = COLOR
    # plt.clf()

    pdf_path = os.path.join(data_path, 'raw_imgs')
    mkdir(pdf_path)
    for constraint_type, fea_history in histories.items():
        unit = '[kN/cm2]' if constraint_type == 'stress' else '[cm]'
        scale = 1.0 if constraint_type == 'stress' else 100.0

        fig, ax = plt.subplots(figsize=(10, 3.5), layout='constrained')
        # figsize=(10, 2.5)

        fea_history = scale * np.array(fea_history)
        sorted_ids = np.argsort(-fea_history)
        max_id = sorted_ids[0]
        max_fea = fea_history[max_id]

        # * performance line
        # perf_color = '#4A3981' if constraint_type == 'displacement' else '#55B679'
        perf_color = 'magenta' if constraint_type == 'displacement' else '#00B0F0'
        ax.plot(range(len(seq)), fea_history, c=perf_color, linewidth = 6)

        # * tolerance line
        # if constraint_type == given_constraint_type:
        #     tol = float(test_folder_name.split('_')[0].split('=')[-1])
        #     ax.plot(range(len(seq)), [tol for _ in seq], c='black', linestyle='--', linewidth = 3)
        #### ax.plot(range(len(seq)), [max_fea for step in seq], c='black', linestyle='-', linewidth = 0.5)

        ax.set_xlabel('assembly steps')
        ax.set_ylabel('maximal {} {}'.format(constraint_type, unit))

        if constraint_type != given_constraint_type:
            # ax.set_ylim(bottom=0.0, top=0.5)
            ax.set_ylim(bottom=0.0, top=0.2)
            # ax.set_ylim(bottom=0.0, top=max_fea + 0.5)
        else:
            # ax.set_ylim(bottom=0.0, top=5.0)
            ax.set_ylim(bottom=0.0, top=1.8)
            # pass

        plt.annotate('{:.2f}'.format(max_fea), xy=(max_id, max_fea), xytext=(max_id, max_fea+0.05), xycoords='data',
                     arrowprops=dict(arrowstyle="->", lw=0.5), fontsize=28,
                     )

        # extraticks = [max_fea]
        # plt.yticks(list(plt.yticks()[0]) + extraticks);

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        fig.savefig(os.path.join(pdf_path, 'design{}_seq{}_{}.png'.format(design_id, seq_id, constraint_type)), dpi=200, transparent=True)
        fig.savefig(os.path.join(pdf_path, 'design{}_seq{}_{}.pdf'.format(design_id, seq_id, constraint_type)), dpi=200)

    return True