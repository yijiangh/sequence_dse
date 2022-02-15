function get_deformation(model, feb, existing_e_ids)
	#    # clone the model and its list of elements to avoid side effects
    #    model = model_in.Clone()
    #    # clone its elements to avoid side effects
    #    model.cloneElements()
    #    # clone the feb-model to avoid side effects
   	# model.deepCloneFEModel()
    
    for e in model.elems
        e.set_is_active(model, e.ind in existing_e_ids)
    end
    
    deform = feb.Deform(model.febmodel)
    response = feb.Response(deform)
    
    try
        # calculate the displacements
        response.updateNodalDisplacements();
        # calculate the member forces
        response.updateMemberForces();
    catch e
        throw("$e, The stiffness matrix of the system is singular.")
    end

    # # if something changed inform the feb-model about it (otherwise it won't recalculate)
    model.febmodel.touch()
    
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
    
   #  this guards the objects from being freed prematurely
   # GC.KeepAlive(deform)
   # GC.KeepAlive(response)
    
	return response.maxDisplacement()
end