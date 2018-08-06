import os
import unittest

from aether.diagnostics import set_diagnostics_on
from aether.geometry import define_node_geometry_2d, define_elem_geometry_2d
from aether.arrays import check_node_xyz_2d, check_node_dxyz_2d,check_elem_2d_version, check_nodes_in_elem_2d

# Look to see if the 'TEST_RESOURCES_DIR' is set otherwise fallback to a path
# relative to this file.
if 'TEST_RESOURCES_DIR' in os.environ:
    resources_dir = os.environ['TEST_RESOURCES_DIR']
else:
    here = os.path.abspath(os.path.dirname(__file__))
    resources_dir = os.path.join(here, 'resources')


class GeometryTestCase(unittest.TestCase):

    def test_read_square(self):
        set_diagnostics_on(False)
        define_node_geometry_2d(os.path.join(resources_dir, 'square.ipnode'))
#       value = check_node_dxyz_2d(2,2, 1)
#       self.assertEqual(-1, value)
#       value = check_node_xyz_2d(1, 1)
#       self.assertEqual(100, value)

#       define_elem_geometry_2d(os.path.join(resources_dir, 'square.ipelem'))

        define_elem_geometry_2d(os.path.join(resources_dir, 'square.ipelem') ,'unit')
        value = check_nodes_in_elem_2d(1, 1)
        self.assertEqual(52, value)

        value = check_elem_2d_version(1, 1)
        self.assertEqual(1, value)


if __name__ == '__main__':
    unittest.main()
