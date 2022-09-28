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
from Karamba.Results import BeamResultantForces, BeamForces
from Karamba.Results import Utilization_Beam, UtilizationResults_Beam
import feb # Karamba's C++ library (undocumented in the API)
#import System.GC as GC
from System.Collections.Generic import List
from System import Guid

##############################

def test_stiffness(model_in, existing_e_ids, lc=0, verbose=False):
    if not existing_e_ids:
        return True
    
    # clone the model and its list of elements to avoid side effects
    model = model_in.Clone()
    # clone its elements to avoid side effects
    model.cloneElements()
    # clone the feb-model to avoid side effects
    model.deepCloneFEModel()
    
    for e in model.elems:
        e.set_is_active(model, e.ind in existing_e_ids)
    
    deform = feb.Deform(model.febmodel)
    response = feb.Response(deform)
    
    try:
        # calculate the displacements
        response.updateNodalDisplacements();
        # calculate the member forces
        response.updateMemberForces();
    except:
        raise ValueError("The stiffness matrix of the system is singular.")
        
    # if something changed inform the feb-model about it (otherwise it won't recalculate)
    model.febmodel.touch()
    
    if not use_stress_constraint:
#        displ_norms = []
        max_disp = -10e-5
        for i in range(len(model.nodes)):
            node_x = model.febmodel.state(0).nodeState(i).disp_trans_global().x()
            node_y = model.febmodel.state(0).nodeState(i).disp_trans_global().y()
            node_z = model.febmodel.state(0).nodeState(i).disp_trans_global().z()
            n_disp = math.sqrt(node_x**2 + node_y**2 + node_z**2)
            if n_disp > max_disp:
                max_disp = n_disp
        return max_disp # m
    else:
        lc = 0
        # no extra identifier allowed for now
        elem_id = ''
        assert len(model.elementsByID(List[str]([elem_id]))) == len(model.elems)
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
        
        max_disp = clr.Reference[List[float]]()
        out_g = clr.Reference[List[float]]()
        out_comp = clr.Reference[List[float]]()
        message = clr.Reference[str]()
        k3d = KarambaCommon.Toolkit();
        model = k3d.Algorithms.AnalyzeThI(model, 
            max_disp, out_g, out_comp, message);
        
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
        return max_sigma
        
#        # https://www.karamba3d.com/help/2-2-0/html/0c37d3de-2cf9-00b7-8b7a-a1ab8734227b.htm
#        lcind = 0
#        out_N = clr.Reference[List[List[float]]]()
#        out_V = clr.Reference[List[List[float]]]()
#        out_M = clr.Reference[List[List[float]]]()
#        
#        # maximum distance between points where section forces are determined.
#        # maximum number of points where section forces are determined. Can be overruled by max_res_dist.
#        BeamResultantForces.solve(model, List[str]([elem_id]),
#            str(lcind), 100, 1, out_N, out_V, out_M)
#        N = out_N.Value # Largest absolute normal force along beam. List structure: load-case/element.
#        V = out_V.Value # Largest resultant shear force along beam
#        M = out_M.Value # Largest resultant bending moment along beam
#        
#        max_sigma = -10e-5
#        for i in range(len(model.elems)):
#            crosec = model.elems[i].crosec
#            mat = crosec.material
#            # round crosec
#            assert abs(abs(crosec.Wely_z_pos) - abs(crosec.Wely_z_neg)) < 1e-10
#            e_max_sigma = abs(N[lcind][i]) / crosec.A + abs(M[lcind][i]) / crosec.Wely_z_pos;
#            if max_sigma < e_max_sigma:
#                max_sigma = e_max_sigma
##            print 'E{}: {:.2f} | N {:.2f}, M {:.2f}, Welyz {}'.format(i, max_sigma, 
##                N[lcind][i], M[lcind][i], crosec.Wely_z_pos)
#        return max_sigma # kN/m2
         
#    energy_visitor = feb.EnergyVisitor(model.febmodel, model.febmodel.state(0), lc);
#    energy_visitor.visit(model.febmodel);
#    elastic_E = energy_visitor.elasticEnergy()
#    if verbose:
#        print '---'
##        print existing_e_ids
#        print 'nKinematicMode: ', deform.nKinematicModes()
#        print 'maxDisp (m): ', response.maxDisplacement()
#        print 'compliance: ', deform.compliance(lc)
#        print 'elastic E (kN-m): ', elastic_E
#        print '---'

#    axial_E = energy_visitor.axialEnergy()
#    bending_E = energy_visitor.bendingEnergy()
#    print 'axial + bending - elastic', abs(axial_E + bending_E - elastic_E)
    
    # this guards the objects from being freed prematurely
#    GC.KeepAlive(deform)
#    GC.KeepAlive(response)
    
#    return response.maxDisplacement()
    
    

def compute_path_cost(plan):
    return max([abs(max_disp) for _, max_disp in plan])

##########################
test_stiffness_fn = test_stiffness
compute_path_cost_fn = compute_path_cost