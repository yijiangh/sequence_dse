function get_deformation_from_karamba(model, feb, existing_e_ids)
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
   
	return response.maxDisplacement()
end