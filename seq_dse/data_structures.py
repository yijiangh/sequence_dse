from collections import namedtuple

SeqElement = namedtuple('Element', ['name', 'fem_data', 'planning_data'])

# used for FEM
FEMData = namedtuple('FEMData', ['axis_points']) # , 'radius'
# used for Pybullet planning
PlanningData = namedtuple('PlanningData', ['start_pose', 'goal_pose', 'body', 'grasp_offsets'])
