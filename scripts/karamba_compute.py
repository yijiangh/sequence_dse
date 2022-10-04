import math
##############################
import os
import clr
# https://ironpython.net/documentation/dotnet/dotnet.html

clr.AddReferenceToFileAndPath(os.path.join(karamba_plugin_path, "Karamba.gha"))
clr.AddReferenceToFileAndPath(os.path.join(karamba_plugin_path, "KarambaCommon.dll"))

import KarambaCommon
import Karamba.Models.Model as Model
import Karamba.GHopper.Models.GH_Model as GH_Model
import Karamba.Elements.ModelTruss as Truss
import Karamba.Elements.ModelBeam as Beam
import Karamba.Materials.FemMaterial_Isotrop as FemMaterial
from Karamba.Geometry import Vector3
#from Karamba.Results import BeamResultantForces, BeamForces
from Karamba.Results import Utilization_Beam, UtilizationResults_Beam
from Karamba.Results import NodalDisp

#import feb # Karamba's C++ library (undocumented in the API)
#import System.GC as GC
from System.Collections.Generic import List
from System import Guid

##############################

def compute_path_cost(plan):
    return max([abs(max_disp if objective=='displacement' else max_stress) \
        for _, (max_disp, max_stress) in plan])

def test_stiffness(model_in, existing_e_ids, lc='in_construction_load', verbose=False):
    if not existing_e_ids:
        return 0, 0
    
    # clone the model and its list of elements to avoid side effects
    model = model_in.Clone()
    # clone its elements to avoid side effects
    model.cloneElements()
    # clone the feb-model to avoid side effects
    model.deepCloneFEModel()
    
    for e in model.elems:
        e.set_is_active(model, e.ind in existing_e_ids)
    
    # * Computing max stress
    # no extra identifier allowed for now
#    elem_id = ''
#    assert len(model.elementsByID(List[str]([elem_id]))) == len(model.elems)
    # number of sampling points along beams
    samplingPointsCount = 5;
    # elastic or plastic cross section design?
    isElasticDesign = True;
    # material safety factor is buckling does not govern the design
    gammaM0 = 1.0;
    # material safety factor is buckling does govern the design
    gammaM1 = 1.1;
    # does buckling involve sideways sway?
    swayFrame = False;
    # shall the output contain calculation details?
    withDetails = False;
    
    max_disp2 = clr.Reference[List[float]]()
    out_g = clr.Reference[List[float]]()
    out_comp = clr.Reference[List[float]]()
    message = clr.Reference[str]()
    k3d = KarambaCommon.Toolkit();
    model = k3d.Algorithms.AnalyzeThI(model, 
        max_disp2, out_g, out_comp, message)
#    print max_disp2
    
    # resulting beam utilizations
    modelUtil = clr.Reference[UtilizationResults_Beam]()
    # resulting message
    message = clr.Reference[str]()
    
    Utilization_Beam.solve(
      model,
      List[str](['']),
      List[Guid]([]),
      str(lc),
      samplingPointsCount,
      isElasticDesign,
      gammaM0,
      gammaM1,
      swayFrame,
      withDetails,
      modelUtil,
      message);
    
    max_sigma = -10e-5
    for i in existing_e_ids: # range(len(model.elems)):
        esigmax = max([abs(modelUtil.sig_max[i]), abs(modelUtil.sig_min[i])])
        if max_sigma < esigmax:
            max_sigma = esigmax
    
    # * Computing max displacement
    disp_trans = clr.Reference[List[List[Vector3]]]()
    disp_rots = clr.Reference[List[List[Vector3]]]()
    node_ids = List[int]([nd.ind for nd in model.nodes])
    NodalDisp.solve(model, lc, node_ids, disp_trans, disp_rots)
    assert len(disp_trans.Value) == 1
    
    max_disp = -10e-5
    for i in range(len(model.nodes)):
        disp_vec = disp_trans.Value[0][i]
        node_x, node_y, node_z = disp_vec.X, disp_vec.Y, disp_vec.Z
        n_disp = math.sqrt(node_x**2 + node_y**2 + node_z**2)
        if n_disp > max_disp:
            max_disp = n_disp
#    print 'max disp {}'.format(max_disp)
    return max_disp, max_sigma
    
##########################
test_stiffness_fn = test_stiffness
compute_path_cost_fn = compute_path_cost