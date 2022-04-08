import numpy as np

def seven_bar_truss_conmech_model(X):
    """
    # O------C-------S
    # \    / | \    /
    #  \  /  |  \  /
    #   \/   |   \/
    #   X--------X'
    # O = (0,0), X = (x,y) the variables
    # norm(O,S) = 2*h_span
    # P is a vertical point load applied at C
    # X' node is mirrored
    """
    # Point load in the middle, kN
    P = -100.0
    # half-span, meter
    h_span = 2.0

    points = np.array([[0.,0.],
                       [h_span,0.],
                       [2*h_span,0.],
                       [X[0], X[1]],
                       [2*h_span-X[0], X[1]]])