import os
from typing import List

from compas.datastructures import Mesh
from compas.geometry import Frame

from incon_building_plan.building_step_object import BuildingStepObject
from incon_building_plan.building_step import (
    ObjectBuildingStep,
    InstructionBuildingStep,
    TagBuildingStep
)
from incon_building_plan.incon_building_plan import InconBuildingPlan
from incon_building_plan.instructions.dialog_2d_instruction import Dialog2dInstruction
from incon_building_plan.instructions.linear_dimension_instruction import (
    LinearDimensionInstruction,
)
from incon_building_plan.instructions.model_instruction import (
    ModelInstruction,
)
from incon_building_plan.actions.clickable_dialog_action import (
    ClickableDialogAction,
)
from incon_building_plan.actions.clickable_video_action import (
    ClickableVideoAction,
)

from .utils import mkdir

def seq_to_incon_instruction(building_plan_id: str, meshes: List[Mesh], mesh_frames: List[Frame], assist_meshes: List[Mesh], 
        assist_mesh_frames: List[Frame], _plan_save_path: str, tag_frames: List[Frame], attached_seq_info: List[str]):
    building_plan = InconBuildingPlan(building_plan_id=building_plan_id)
    plan_save_path = os.path.join(_plan_save_path, building_plan_id)
    mkdir(plan_save_path)
    mesh_save_path = os.path.join(plan_save_path, 'meshes')
    mkdir(mesh_save_path)

    # * add AprilTag building steps
    for i, frame in enumerate(tag_frames):
        point = frame.point
        quat = frame.quaternion 
        tag_step = TagBuildingStep(
            id=f'tag_{i}',
            tag_id=i,
            tag_size=0.094, # if the tag is printed unscaled, directly from the PDF
            pos_x=point.x,
            pos_y=point.y,
            pos_z=point.z,
            quat_x=quat.x,
            quat_y=quat.y,
            quat_z=quat.z,
            quat_w=quat.w,
        )
        building_plan.add_building_step(tag_step)

    # * add registeration-asist element building steps
    for i, (mesh, mesh_frame) in enumerate(zip(assist_meshes, assist_mesh_frames)):
        point = mesh_frame.point
        quat = mesh_frame.quaternion 
        # * save mesh to a path
        mesh_file_name = f"assist_element_{i}.obj"
        mesh.to_obj(os.path.join(mesh_save_path,mesh_file_name))

        # * create one building step
        building_step = ObjectBuildingStep(
            build_instructions=[],
            color_rgb=[0.0, 0.0, 1.0],
            dependencies=[],
            id=f"assist_element_{i}",
            is_already_built=True,
            object_type=mesh_file_name,
            pos_x=point.x,
            pos_y=point.y,
            pos_z=point.z,
            quat_x=quat.x,
            quat_y=quat.y,
            quat_z=quat.z,
            quat_w=quat.w,
        )
        building_plan.add_building_step(building_step)

    # * add real element building steps
    for i, (mesh, mesh_frame, info) in enumerate(zip(meshes, mesh_frames, attached_seq_info)):
        # * instructions
        dialog_2d_instruction = Dialog2dInstruction(
            id=f"element_{i}_instruction_2d_diag",
            title=f'Seq#{i}',
            text=info,
        )

        point = mesh_frame.point
        quat = mesh_frame.quaternion 

        # * save mesh to a path
        mesh_file_name = f"element_{i}.obj"
        mesh.to_obj(os.path.join(mesh_save_path,mesh_file_name))

        # * create one building step
        building_step = ObjectBuildingStep(
            build_instructions=[],
            color_rgb=[1.0, 1.0, 1.0],
            dependencies=[],
            id=f"element_{i}",
            is_already_built=False,
            object_type=mesh_file_name,
            pos_x=point.x,
            pos_y=point.y,
            pos_z=point.z,
            quat_x=quat.x,
            quat_y=quat.y,
            quat_z=quat.z,
            quat_w=quat.w,
        )
        building_step.add_dialog_instruction(dialog_2d_instruction)
        building_plan.add_building_step(building_step)

    building_plan_path = os.path.join(plan_save_path, building_plan_id + ".json")
    print("Writing building plan to: " + building_plan_path)
    building_plan.write_to_file(path=building_plan_path)