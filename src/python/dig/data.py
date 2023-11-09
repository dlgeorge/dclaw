#!/usr/bin/env python

"""

Classes representing parameters for GeoClaw runs

:Classes:

 - GeoClawData
 - RefinementData
 - TopographyData
 - FixedGridData
 - FGmaxData
 - DTopoData
 - QinitData
 - SurgeData
 - MultilayerData
 - FrictionData 
 - GridData1D
 - BoussData1D

:Constants:

 - Rearth - Radius of earth in meters
 - DEG2RAD factor to convert degrees to radians
 - RAD2DEG factor to convert radians to degrees
 - LAT2METER factor to convert degrees in latitude to meters
"""

from __future__ import absolute_import
from __future__ import print_function
import os
import numpy
import clawpack.clawutil.data
import warnings


# Radius of earth in meters.
# For consistency, should always use this value when needed, e.g.
# in setrun.py or topotools:
Rearth = 6367.5e3  # average of polar and equatorial radii

DEG2RAD = numpy.pi / 180.0
RAD2DEG = 180.0 / numpy.pi
LAT2METER = Rearth * DEG2RAD

class GeoClawData(clawpack.clawutil.data.ClawData):
    r"""
    Object containing the basic .

    Note that this data object will write out multiple files.
    """
    def __init__(self):
        super(GeoClawData,self).__init__()

        # GeoClaw physics parameters
        self.add_attribute('gravity',9.8)
        self.add_attribute('rho', 1025.0)  # Density of water kg/m^3
        self.add_attribute('rho_air',1.15) # Density of air kg/m^3
        self.add_attribute('ambient_pressure', 101.3e3) # Nominal atmos pressure
        self.add_attribute('earth_radius',Rearth)
        self.add_attribute('coordinate_system',1)
        self.add_attribute('sphere_source',0)  # should set to 1 by default?
        self.add_attribute('coriolis_forcing',True)
        self.add_attribute('theta_0',45.0)
        self.add_attribute('friction_forcing',True)
        self.add_attribute('manning_coefficient',[0.025])
        self.add_attribute('manning_break',[])

        # GeoClaw algorithm parameters
        self.add_attribute('dry_tolerance',1e-3)
        self.add_attribute('friction_depth',1.0e6)
        self.add_attribute('sea_level',0.0)


    def write(self,data_source='setrun.py', out_file='geoclaw.data'):

        self.open_data_file(out_file, data_source)

        self.self.data_write('gravity',
                               description="(gravitational acceleration m/s^2)")
        self.self.data_write('rho', description="(Density of water kg/m^3)")
        self.self.data_write('rho_air',description="(Density of air kg/m^3)")
        self.self.data_write('ambient_pressure',
                                description="(Nominal atmospheric pressure Pa)")
        self.self.data_write('earth_radius', description="(Radius of the earth m)")
        self.self.data_write('coordinate_system',
                        description="(1=meters, 2=lon-lat)")
        self.self.data_write('sphere_source',
                        description="(0=none, 1=only in mass eqn, 2=all)")
        self.self.data_write('sea_level')
        self.self.data_write()

        # Forcing terms
        self.self.data_write('coriolis_forcing')
        if self.coordinate_system == 1 and self.coriolis_forcing:
            self.self.data_write('theta_0')
        self.self.data_write('friction_forcing')
        if self.friction_forcing:
            if type(self.manning_coefficient) in [int,float]:
                self.manning_coefficient = [self.manning_coefficient]
            num_manning = len(self.manning_coefficient)
            if len(self.manning_break) != num_manning - 1:
                raise IOError("***manning_break array has wrong length")
            self.self.data_write(value=num_manning,alt_name='num_manning')
            self.self.data_write('manning_coefficient')
            self.self.data_write('manning_break')
            self.self.data_write('friction_depth')

        self.self.data_write()

        self.self.data_write('dry_tolerance')

        self.close_data_file()



class RefinementData(clawpack.clawutil.data.ClawData):

    def __init__(self):

        super(RefinementData,self).__init__()

        # Refinement controls
        self.add_attribute('wave_tolerance',1.0e-1)
        self.add_attribute('speed_tolerance',[1.0e12]*6)
        self.add_attribute('deep_depth',None)      # deprecated
        self.add_attribute('max_level_deep',None)  # deprecated
        self.add_attribute('variable_dt_refinement_ratios',False)


    def write(self,data_source='setrun.py', out_file='refinement.data'):
        # Refinement controls
        self.open_data_file(out_file, data_source)
        self.self.data_write('wave_tolerance')

        # check if user set deprecated parameters:
        if self.deep_depth is not None:
            w = '\n  *** WARNING: deep_depth parameter ignored as of v5.8.0'
            warnings.warn(w, UserWarning)
        if self.max_level_deep is not None:
            w = '\n  *** WARNING: max_level_deep parameter ignored as of v5.8.0'
            warnings.warn(w, UserWarning)

        if not isinstance(self.speed_tolerance,list):
            self.speed_tolerance = [self.speed_tolerance]
        self.self.data_write('speed_tolerance')
        self.self.data_write()
        self.self.data_write('variable_dt_refinement_ratios',
                        description="(Set dt refinement ratios automatically)")
        self.close_data_file()



class TopographyData(clawpack.clawutil.data.ClawData):

    def __init__(self):

        super(TopographyData,self).__init__()

        # Topography data
        self.add_attribute('topo_missing',99999.)
        self.add_attribute('test_topography',0)
        self.add_attribute('topofiles',[])

        # Jump discontinuity
        self.add_attribute('topo_location',-50e3)
        self.add_attribute('topo_left',-4000.0)
        self.add_attribute('topo_right',-200.0)
        self.add_attribute('topo_angle',0.0)

        # Simple oceanic shelf
        self.add_attribute('x0',350e3)
        self.add_attribute('x1',450e3)
        self.add_attribute('x2',480e3)
        self.add_attribute('basin_depth',-3000.0)
        self.add_attribute('shelf_depth',-100.0)
        self.add_attribute('beach_slope',0.008)


    def write(self,data_source='setrun.py', out_file='topo.data'):

        self.open_data_file(out_file, data_source)
        self.self.data_write(name='topo_missing',
                        description='replace no_data_value in topofile')
        self.self.data_write(name='test_topography',description='(Type topography specification)')
        if self.test_topography == 0:
            ntopofiles = len(self.topofiles)
            self.self.data_write(value=ntopofiles,alt_name='ntopofiles')
            for tfile in self.topofiles:

                if len(tfile) == 6:
                    w = '\n  *** WARNING: topofile specs changed in v5.8.0 -- ' + \
                          'Flag level info now ignored'
                    warnings.warn(w, UserWarning)
                    tfile = [tfile[0], tfile[-1]] # drop minlevel,maxlevel,t1,t2
                elif len(tfile) == 2:
                    pass  # now expect only topo_type, filename
                else:
                    raise ValueError('Unexpected len(tfile) = %i' % len(tfile))

                # if path is relative in setrun, assume it's relative to the
                # same directory that out_file comes from
                fname = os.path.abspath(os.path.join(os.path.dirname(out_file),tfile[-1]))
                self._out_file.write("\n'%s' \n " % fname)
                self._out_file.write("%3i   # topo_type\n" % tfile[0])
        elif self.test_topography == 1:
            self.self.data_write(name='topo_location',description='(Bathymetry jump location)')
            self.self.data_write(name='topo_left',description='(Depth to left of bathy_location)')
            self.self.data_write(name='topo_right',description='(Depth to right of bathy_location)')
        elif self.test_topography == 2 or self.test_topography == 3:
            self.self.data_write(name='x0',description='(Location of basin end)')
            self.self.data_write(name='x1',description='(Location of shelf slope end)')
            self.self.data_write(name='x2',description='(Location of beach slope)')
            self.self.data_write(name='basin_depth',description='(Depth of basin)')
            self.self.data_write(name='shelf_depth',description='(Depth of shelf)')
            self.self.data_write(name='beach_slope',description='(Slope of beach)')
        else:
            raise NotImplementedError("Test topography type %s has not been"
                                        " implemented." % self.test_topography)

        self.close_data_file()


class FixedGridData(clawpack.clawutil.data.ClawData):

    """
    Deprecated, starting in 5.9.0 use FGoutData instead.
    """

    def __init__(self):

        super(FixedGridData,self).__init__()

        # Fixed Grids
        self.add_attribute('fixedgrids',[])


    def write(self,data_source='setrun.py', out_file='fixed_grids.data'):
        # Fixed grid settings
        msg = 'rundata.fixed_grid_data is deprecated starting in v5.9.0,' \
            + ' use rundata.fgout_data instead'
        #warnings.warn(msg)
        if len(self.fixedgrids) > 0:
            raise AttributeError(msg)


class FGoutData(clawpack.clawutil.data.ClawData):

    def __init__(self):

        super(FGoutData,self).__init__()

        # File name for fgout points and parameters:
        self.add_attribute('fgout_grids',[])


    def write(self,data_source='setrun.py', out_file='fgout_grids.data'):
        self.open_data_file(out_file, data_source)
        num_fgout_grids = len(self.fgout_grids)
        self.self.data_write(value=num_fgout_grids,alt_name='num_fgout_grids')
        self.self.data_write()

        fgno_unset = 0  # to use if fg.fgno not set by user
        fgno_list = []  # to check for uniqueness of fgno's

        for fg in self.fgout_grids:
            # if path is relative in setrun, assume it's relative to the
            # same directory that out_file comes from

            if fg.fgno is None:
                # not set by user in setrun
                fgno_unset += 1
                fg.fgno = fgno_unset

            if fg.fgno in fgno_list:
                msg = 'Trying to set fgout grid number to fgno = %i' % fg.fgno \
                      + '\n             but this fgno was already used' \
                      + '\n             Set unique fgno for each fgout grid'
                raise ValueError(msg)

            fgno_list.append(fg.fgno)
            fg.write_to_fgout_data(self._out_file)
        self.close_data_file()

class FGmaxData(clawpack.clawutil.data.ClawData):

    def __init__(self):

        super(FGmaxData,self).__init__()

        # File name for fgmax points and parameters:
        self.add_attribute('fgmax_files',[])
        self.add_attribute('num_fgmax_val',1)
        self.add_attribute('fgmax_grids',[])


    def write(self,data_source='setrun.py', out_file='fgmax_grids.data'):
        if len(self.fgmax_files) > 0:
            msg = '*** fgmax_files has been deprecated, ' \
                  + 'use fgmax_grids instead.'
            raise ValueError(msg)

        # new style:
        self.open_data_file(out_file, data_source)
        num_fgmax_val = self.num_fgmax_val
        if num_fgmax_val not in [1,2,5]:
            raise NotImplementedError(
                   "Expecting num_fgmax_val in [1,2,5], got %s" % num_fgmax_val)
        self.self.data_write(value=num_fgmax_val, alt_name='num_fgmax_val')
        num_fgmax_grids = len(self.fgmax_grids)
        self.self.data_write(value=num_fgmax_grids, alt_name='num_fgmax_grids')
        self.self.data_write()

        fgno_unset = 0  # to use if fg.fgno not set by user
        fgno_list = []  # to check for uniqueness of fgno's

        for fg in self.fgmax_grids:
            # if path is relative in setrun, assume it's relative to the
            # same directory that out_file comes from
            if fg.xy_fname is not None:
                fg.xy_fname = os.path.abspath(os.path.join(\
                              os.path.dirname(out_file),fg.xy_fname))

            if fg.fgno is None:
                # not set by user in setrun
                fgno_unset += 1
                fg.fgno = fgno_unset

            if fg.fgno in fgno_list:
                msg = 'Trying to set fgmax grid number to fgno = %i' % fg.fgno \
                      + '\n             but this fgno was already used' \
                      + '\n             Set unique fgno for each fgmax grid'
                raise ValueError(msg)

            fgno_list.append(fg.fgno)
            fg.write_to_fgmax_data(self._out_file)
        self.close_data_file()


    def read(self, path="fgmax_grids.data", force=False):
        r"""Read a FGMax data file."""

        super(FGmaxData, self).read(path, force=force)

        # Look for basic parameters
        fig_numbers = []
        with open(os.path.abspath(path), 'r') as data_file:
            # Forward to first parameter
            for line in data_file:
                # Regular parameter setting
                if "=:" in line:
                    value, tail = line.split("=:")
                    varname = tail.split()[0]

                    if varname == "num_fgmax_val":
                        self.num_fgmax_val = int(value)
                    elif varname == "num_fgmax_grids":
                        num_fgmax_grids = int(value)
                
                # Contains a fixed grid number
                elif "# fgno" in line:
                    value, tail = line.split("#")
                    fig_numbers.append(int(value))

        if len(fig_numbers) != num_fgmax_grids:
            raise ValueError("Number of FGMaxGrid numbers found does not ", 
                             "equal the number of grids recorded.")
        
        # Read each fgmax grid
        import clawpack.geoclaw.fgmax_tools
        for (i, grid_num) in enumerate(fig_numbers):
            new_fgmax_grid = clawpack.geoclaw.fgmax_tools.FGmaxGrid()
            new_fgmax_grid.read_fgmax_grids_data(grid_num, data_file=path)
            self.fgmax_grids.append(new_fgmax_grid)



class DTopoData(clawpack.clawutil.data.ClawData):

    def __init__(self):

        super(DTopoData,self).__init__()

        # Moving topograhpy
        self.add_attribute('dtopofiles',[])
        self.add_attribute('dt_max_dtopo', 1.e99)

    def write(self, data_source='setrun.py', out_file='dtopo.data'):

        # Moving topography settings
        self.open_data_file(out_file, data_source)
        mdtopofiles = len(self.dtopofiles)
        self.self.data_write(value=mdtopofiles,alt_name='mdtopofiles')
        self.self.data_write()
        for tfile in self.dtopofiles:

            if len(tfile) == 4:
                w = '\n  *** WARNING: dtopofile specs changed in v5.8.0 -- ' + \
                      'Flag level info now ignored'
                warnings.warn(w, UserWarning)
                tfile = [tfile[0], tfile[-1]]  # drop minlevel,maxlevel
            elif len(tfile) == 2:
                pass   # now expect only dtopo_type, filename
            else:
                raise ValueError('Unexpected len(tfile) = %i' % len(tfile))

            # if path is relative in setrun, assume it's relative to the
            # same directory that out_file comes from
            fname = os.path.abspath(os.path.join(os.path.dirname(out_file),tfile[-1]))
            self._out_file.write("\n'%s' \n" % fname)
            self._out_file.write("%3i   # dtopo_type\n" % tfile[0])
        self.self.data_write()
        self.self.data_write(value=self.dt_max_dtopo,alt_name='dt_max_dtopo')
        self.close_data_file()


    def read(self, path="dtopo.data", force=False):
        r"""Read a dtopography data file."""

        print(self.dtopofiles)

        with open(os.path.abspath(path), 'r') as data_file:

            file_name = None

            # Forward to first parameter
            for line in data_file:

                # Regular parameter setting
                if "=:" in line:
                    value, tail = line.split("=:")
                    varname = tail.split()[0]

                    if varname == "mdtopofiles":
                        num_dtopo_files = int(value)
                    elif varname == "dt_max_dtopo":
                        self.dt_max_dtopo = float(value)

                # Assume this is the second line of a record, complete and add
                # to dtopofiles list
                elif file_name is not None:
                    base_values = [int(value) for value in line.split()]
                    base_values.append(file_name)
                    self.dtopofiles.append(base_values)
                    file_name = None

                # Non-empty line, assume start of dtopo_file record
                elif line[0] == "'":
                    file_name = line.strip()[1:-1]

        # Check to make sure we have all the dtopofiles
        if len(self.dtopofiles) != num_dtopo_files:
            raise IOError("The number of dtopo files specified does not equal ",
                          "the number found.")



class ForceDry(clawpack.clawutil.data.ClawData):

    def __init__(self):
        r"""
        A single force_dry array and associated data
        """

        super(ForceDry,self).__init__()
        self.add_attribute('tend',None)
        self.add_attribute('fname','')


class QinitData(clawpack.clawutil.data.ClawData):

    def __init__(self):

        super(QinitData,self).__init__()

        # Qinit data
        self.add_attribute('qinit_type',0)
        self.add_attribute('qinitfiles',[])
        self.add_attribute('variable_eta_init',False)
        self.add_attribute('force_dry_list',[])
        self.add_attribute('num_force_dry',0)

    def write(self,data_source='setrun.py', out_file='qinit.data'):

        # Initial perturbation
        self.open_data_file(out_file, data_source)
        self.self.data_write('qinit_type')

        # Perturbation requested
        if self.qinit_type == 0:
            pass
        else:
            # Check to see if each qinit file is present and then write the data
            for tfile in self.qinitfiles:

                if len(tfile) == 3:
                    w = '\n  *** WARNING: qinit specs changed in v5.8.0 -- ' + \
                          'Flag level info now ignored'
                    warnings.warn(w, UserWarning)
                    tfile = [tfile[-1]]  # drop minlevel,maxlevel
                elif len(tfile) == 1:
                    pass  # now expect only filename
                else:
                    raise ValueError('Unexpected len(tfile) = %i' % len(tfile))

                # if path is relative in setrun, assume it's relative to the
                # same directory that out_file comes from
                fname = os.path.abspath(os.path.join(os.path.dirname(out_file),tfile[-1]))
                self._out_file.write("\n'%s' \n" % fname)
        # else:
        #     raise ValueError("Invalid qinit_type parameter %s." % self.qinit_type)


        self.self.data_write('variable_eta_init')

        self.num_force_dry = len(self.force_dry_list)
        self.self.data_write('num_force_dry')

        for force_dry in self.force_dry_list:

            # if path is relative in setrun, assume it's relative to the
            # same directory that out_file comes from
            fname = os.path.abspath(os.path.join(os.path.dirname(out_file),\
                    force_dry.fname))
            self._out_file.write("\n'%s' \n" % fname)
            self._out_file.write("%.3f \n" % force_dry.tend)


        self.close_data_file()


# Storm data
class SurgeData(clawpack.clawutil.data.ClawData):
    r"""Data object describing storm surge related parameters"""

    # Provide some mapping between model names and integers
    storm_spec_dict_mapping = {"HWRF":-1,
                               None: 0,
                               'holland80': 1,
                               'holland08': 8,
                               'holland10': 2,
                               'CLE': 3,
                               'SLOSH': 4,
                               'rankine': 5,
                               'modified-rankine': 6,
                               'DeMaria': 7
                              }
    storm_spec_not_implemented = ['CLE']

    def __init__(self):
        super(SurgeData,self).__init__()

        # Source term controls
        self.add_attribute('wind_forcing',False)
        self.add_attribute('drag_law',1)
        self.add_attribute('pressure_forcing',False)

        # Algorithm parameters - Indexing is python based
        self.add_attribute("wind_index", 4)
        self.add_attribute("pressure_index", 6)
        self.add_attribute("display_landfall_time", False)

        # AMR parameters
        self.add_attribute('wind_refine',[20.0,40.0,60.0])
        self.add_attribute('R_refine',[60.0e3,40e3,20e3])

        # Storm parameters
        self.add_attribute('storm_type', None)  # Backwards compatibility
        self.add_attribute('storm_specification_type', 0) # Type of parameterized storm
        self.add_attribute("storm_file", None) # File(s) containing data


    def write(self,out_file='surge.data',data_source="setrun.py"):
        """Write out the data file to the path given"""

        # print "Creating data file %s" % out_file
        self.open_data_file(out_file,data_source)

        self.self.data_write('wind_forcing', description='(Wind source term used)')
        self.self.data_write('drag_law', description='(Type of drag law to use)')
        self.self.data_write('pressure_forcing',
                        description="(Pressure source term used)")
        self.self.data_write()

        self.self.data_write("wind_index", value=self.wind_index + 1,
                        description="(Index into aux array - fortran indexing)")
        self.self.data_write("pressure_index", value=self.pressure_index +  1,
                        description="(Index into aux array - fortran indexing)")
        self.self.data_write("display_landfall_time",
                        description='(Display time relative to landfall)')
        self.self.data_write()

        if isinstance(self.wind_refine, bool):
            if not self.wind_refine:
                self.self.data_write('wind_refine', value=False,
                                description='(Refinement ratios)')
        elif isinstance(self.wind_refine, type(None)):
            self.self.data_write('wind_refine', value=False,
                            description='(Refinement ratios)')
        else:
            self.self.data_write('wind_refine',description='(Refinement ratios)')
        if isinstance(self.R_refine, bool):
            if not self.R_refine:
                self.self.data_write('R_refine', value=False,
                                description='(Refinement ratios)')
        elif isinstance(self.R_refine, type(None)):
            self.self.data_write('R_refine', value=False,
                            description='(Refinement ratios)')
        else:
            self.self.data_write('R_refine', description='(Refinement ratios)')
        self.self.data_write()

        # Storm specification
        if self.storm_type is not None:
            self.storm_specification_type = self.storm_type
        if type(self.storm_specification_type) is not int:
            if self.storm_specification_type in         \
                    self.storm_spec_dict_mapping.keys():
                if self.storm_specification_type in     \
                    self.storm_spec_not_implemented:
                    raise NotImplementedError("%s has not been implemented."
                                %self.storm_specification_type)

                else:
                    self.self.data_write("storm_specification_type",
                                self.storm_spec_dict_mapping[
                                        self.storm_specification_type],
                                description="(Storm specification)")
            else:
                raise ValueError("Unknown storm specification type %s"
                                 % self.storm_specification_type)
        else:
            self.self.data_write("storm_specification_type",
                            description="(Storm specification)")
        self.self.data_write("storm_file", description='(Path to storm data)')

        self.close_data_file()


class FrictionData(clawpack.clawutil.data.ClawData):
    r"""Data class representing complex variable friction"""

    def __init__(self):
        r""""""

        super(FrictionData, self).__init__()

        # Variable friction support
        self.add_attribute('variable_friction', False)

        # Index where the variable friction field is stored (Python indexed)
        self.add_attribute('friction_index', 3)

        # Region support
        self.add_attribute('friction_regions', [])

        # File support
        self.add_attribute('friction_files', [])

    def write(self, out_file='friction.data', data_source='setrun.py'):

        self.open_data_file(out_file, data_source)

        self.self.data_write('variable_friction',
                        description="(method for setting variable friction)")
        self.self.data_write('friction_index', value=self.friction_index + 1,
                        description=("(Index into aux array ",
                                     "- fortran indexing)"))
        self.self.data_write()
        if self.variable_friction:
            # Region based friction
            self.self.data_write(value=len(self.friction_regions),
                            alt_name='num_friction_regions',
                            description="(Friction Regions)")
            self.self.data_write()
            for region in self.friction_regions:
                self.self.data_write(value=region[0], alt_name="lower")
                self.self.data_write(value=region[1], alt_name="upper")
                self.self.data_write(value=region[2], alt_name="depths")
                self.self.data_write(value=region[3],
                                alt_name="manning_coefficients")
                self.self.data_write()

            # File based friction
            self.self.data_write(value=len(self.friction_files),
                            alt_name='num_friction_files')
            for friction_file in self.friction_files:
                # if path is relative in setrun, assume it's relative to the
                # same directory that out_file comes from
                fname = os.path.abspath(os.path.join(os.path.dirname(out_file),friction_file))
                self._out_file.write("'%s' %s\n " % fname)

        self.close_data_file()


class MultilayerData(clawpack.clawutil.data.ClawData):
    r"""
    Multilayer SWE data object
    """

    def __init__(self):
        super(MultilayerData, self).__init__()

        # Physics parameters
        self.add_attribute('num_layers', 1)
        self.add_attribute('rho', [1025.0, 1028.0])
        self.add_attribute('eta', [0.0, -200.0])
        self.add_attribute('wave_tolerance', [1.e-1, 1.e-1])

        # Algorithm parameters
        self.add_attribute('eigen_method', 4)
        self.add_attribute('inundation_method', 2)
        self.add_attribute('check_richardson', True)
        self.add_attribute('richardson_tolerance', 0.95)
        self.add_attribute('layer_index', 8)

        # Need to adjust refinement module for this, dry_limit is in geodata
        self.add_attribute('wave_tolerance', [1e-1, 2e-1])
        self.add_attribute('dry_limit', False)

    def write(self, out_file='multilayer.data', datasource="setrun.py"):

        self.open_data_file(out_file, datasource)

        self.self.data_write('num_layers', description='(Number of layers)')
        self.self.data_write('eta',
                        description='(Initial top surface of each layer)')
        self.self.data_write('wave_tolerance',
                        description=('(Tolerance of surface perturbation per',
                                     ' layer, used for refinement criteria)'))
        self.self.data_write('layer_index', value=self.layer_index + 1,
                        description=("(Index into aux array -",
                                     " fortran indexing)"))
        self.self.data_write(None)
        self.self.data_write('check_richardson',
                        description="(Check Richardson number)")
        self.self.data_write('richardson_tolerance',
                        description='(Tolerance for Richardson number)')
        self.self.data_write('eigen_method',
                        description='(Method for calculating eigenspace)')
        self.self.data_write('inundation_method',
                        description=('(Method for calculating inundation ',
                                     'eigenspace)'))
        self.close_data_file()





#  Gauge data object removed, version from amrclaw works in 1d
#class GaugeData1D(clawpack.clawutil.data.ClawData):


class GridData1D(clawpack.clawutil.data.ClawData):
    r"""
    1D data object for grid info

    """
    def __init__(self):
        super(GridData1D,self).__init__()

        self.add_attribute('grid_type',0)
        self.add_attribute('fname_celledges',None)
        self.add_attribute('monitor_fgmax',False)
        self.add_attribute('monitor_runup',False)
        self.add_attribute('monitor_total_zeta',False)

    def write(self,out_file='grid.data',data_source='setrun.py'):

        self.open_data_file(out_file,data_source)

        self.self.data_write('grid_type')
        if self.grid_type == 2:
            if self.fname_celledges is None:
                self.fname_celledges = 'celledges.txt'
                print('*** grid_type ==2 and fname_celledges not specified,')
                print('*** using celledges.txt')
            # if path is relative in setrun, assume it's relative to the
            # same directory that out_file comes from
            fname = os.path.abspath(os.path.join(os.path.dirname(out_file),
                                    self.fname_celledges))
            self._out_file.write("\n'%s'   =: fname_celledges\n " % fname)

        self._out_file.write("\n%s   =: monitor_fgmax" \
                             % str(self.monitor_fgmax)[0])
        self._out_file.write("\n%s   =: monitor_runup" \
                             % str(self.monitor_runup)[0])
        self._out_file.write("\n%s   =: monitor_total_zeta" \
                             % str(self.monitor_total_zeta)[0])
        self.close_data_file()

    def read(self, path, force=False):
        with open(os.path.abspath(path), 'r') as data_file:
            for line in data_file:
                if "=:" in line:
                    value, tail = line.split("=:")
                    varname = tail.split()[0]
                    if varname == 'grid_type':
                        self.grid_type = int(value)
                    elif varname == 'fname_celledges':
                        self.fname_celledges = value.strip()


class BoussData1D(clawpack.clawutil.data.ClawData):
    r"""
    1D data object for Boussinesq info

    """
    def __init__(self):
        super(BoussData1D,self).__init__()

        self.add_attribute('bouss_equations',2)
        self.add_attribute('bouss_min_depth',20.)

    def write(self,out_file='bouss.data',data_source='setrun.py'):

        self.open_data_file(out_file,data_source)

        self.self.data_write('bouss_equations')
        self.self.data_write('bouss_min_depth')

        self.close_data_file()



# place DIG attributes and description in one place (since this is documents )
# This is the best place to see documentation for D-Claw specific input parameters.
_DIG_ATTRS = {
    "rho_s": "solid grain density (kg/m^3)",
    "rho_f": "pore-fluid density  (kg/m^3)",
    "phi_bed": "basal friction angle (degrees)",
    "theta_input": "slope angle (degrees)",
    "delta": "characteristic grain diameter (m)",
    "kappita": "permeability at m=setdig.m0 (m^2), k0 in G&I eq 2.7 if m0 is 0.6",
    "mu": "viscosity of pore-fluid (Pa-s)",
    "alpha_c": "debris compressibility constant (#)",
    "m_crit": "critical state value of m (#)",
    "c1": "dilation regularization coefficient 1 (#)",
    "m0": "initial solid volume fraction (#)",
    "sigma_0": "baseline stress for definition of compressibility",
    "alpha_seg": "coefficient of segregation velocity profile. When alpha_seg = 0, no segregation occurs",
    "bed_normal": "use of bed normal coordinates (0=false, 1=true). bed_normal = 1 requires theta in aux for slope in one direction",
    "phi_seg_coeff": "adjustment to friction coefficient based on segregation", # not currently used.
    "entrainment": "flag for entrainment, 0 = no entrainment",
    "entrainment_rate": "rate of entrainment parameter 0-1",
    "mom_autostop": "flag for momentum autostop F = no autostop, T = autostop",
    "mom_perc": "percentage of max momentum for autostop, default is 0.05 (5%)",
    "momlevel": "level to do momentum calculation IF mom_autostop==True",
    "src_ftn_num": "number of in-domain sources, if used the file 'sethydrographs.data' is required",
    "fric_offset_val": "start/stop friction offset (degrees). if this value is >0, then hysteretic friction is used (Rocha, Johnson, Gray, 2019)",
    "fric_star_val": "deep friction offset (degrees). only used when fric_offset_val > 0 (Rocha, Johnson, Gray, 2019)",
    "chi_init_val": "initial fraction of species 1, (#). Between 0-1.",
    "kappita_diff": "permeability multiplier for different size species. Only used when alpha_seg>0. kappita is used for species1, kappita*kappita_diff used for species2",
    "outaux": "flag for writing aux to output F = not written, T = written",
    "curvature": "flag for curvature correction 0 = not used, 1 = used",
    "init_ptype": "-1 = zero pressure or user defined files in qinit, 0 = hydrostatic, 1,2 = failure pressure (1=min, 2=avg), 3,4 = rising pressure (3=min, 4=avg)",
    "init_pmax_ratio": "p(init_ptf2)= hydro*init_pmax_ratio: pressure will rise to hydrostatic *init_pmax_ratio",
    "init_ptf": " p(init_ptf) = failure, pressure will rise until t = init_ptf without dilatancy",
    "init_ptf2": "p(init_ptf2)= hydro*init_pmax_ratio, pressure will rise until t = init_ptf2",
}


class DClawInputData(clawpack.clawutil.data.ClawData):
    r"""
    D-Claw data object
    """

    def __init__(self, ndim):
        super(DClawInputData, self).__init__()

        # Set default values:
        self.add_attribute("rho_s", 2700.0)
        self.add_attribute("rho_f", 1000.0)
        self.add_attribute("phi_bed", 40.0)
        self.add_attribute("theta_input")
        self.add_attribute("delta", 0.01)
        self.add_attribute("kappita", 0.0001)
        self.add_attribute("mu", 0.001)
        self.add_attribute("alpha_c", 1.0)
        self.add_attribute("m_crit", 0.62)
        self.add_attribute("c1", 1.0)
        self.add_attribute("m0", 0.52)
        self.add_attribute("sigma_0", 1.0e3)
        self.add_attribute("alpha_seg", 0.0)
        self.add_attribute("bed_normal", 0)
        self.add_attribute("phi_seg_coeff")
        self.add_attribute("entrainment", 0)
        self.add_attribute("entrainment_rate", 0.2)
        self.add_attribute("mom_autostop", False)
        self.add_attribute("mom_perc", 0.05)
        self.add_attribute("src_ftn_num", 0)
        self.add_attribute("fric_offset_val", 0.0)
        self.add_attribute("fric_star_val", 0.0)
        self.add_attribute("chi_init_val", 0.0)
        self.add_attribute("kappita_diff", 1.0)
        self.add_attribute("outaux", False)
        self.add_attribute("curvature", 1)
        self.add_attribute("momlevel", 1)


    def write(self,out_file='setdig.data',data_source='setrun.py'):
        self.open_data_file(out_file,data_source)

        self.self.data_write("rho_s", description=_DIG_ATTRS["rho_s"])
        self.self.data_write("rho_f", description=_DIG_ATTRS["rho_f"])
        self.self.data_write("phi_bed", description=_DIG_ATTRS["phi_bed"])
        self.self.data_write("theta_input", description=_DIG_ATTRS["theta_input"])
        self.self.data_write("delta", description=_DIG_ATTRS["delta"])
        self.self.data_write("kappita", description=_DIG_ATTRS["kappita"])
        self.self.data_write("mu", description=_DIG_ATTRS["mu"])
        self.self.data_write("alpha_c", description=_DIG_ATTRS["alpha_c"])
        self.self.data_write("m_crit", description=_DIG_ATTRS["m_crit"])
        self.self.data_write("c1", description=_DIG_ATTRS["c1"])
        self.self.data_write("m0", description=_DIG_ATTRS["m0"])
        self.self.data_write("sigma_0", description=_DIG_ATTRS["sigma_0"])
        self.self.data_write("alpha_seg", description=_DIG_ATTRS["alpha_seg"])
        self.self.data_write("bed_normal", description=_DIG_ATTRS["bed_normal"])
        self.self.data_write("phi_seg_coeff", description=_DIG_ATTRS["phi_seg_coeff"])
        self.self.data_write("entrainment", description=_DIG_ATTRS["entrainment"])
        self.self.data_write("entrainment_rate", description=_DIG_ATTRS["entrainment_rate"])
        self.self.data_write("mom_autostop", description=_DIG_ATTRS["mom_autostop"])
        self.self.data_write("mom_perc", description=_DIG_ATTRS["mom_perc"])
        self.self.data_write("src_ftn_num", description=_DIG_ATTRS["src_ftn_num"])
        self.self.data_write("fric_offset_val", description=_DIG_ATTRS["fric_offset_val"])
        self.self.data_write("fric_star_val", description=_DIG_ATTRS["fric_star_val"])
        self.self.data_write("chi_init_val", description=_DIG_ATTRS["chi_init_val"])
        self.self.data_write("kappita_diff", description=_DIG_ATTRS["kappita_diff"])
        self.self.data_write("outaux", description=_DIG_ATTRS["outaux"])
        self.self.data_write("curvature", description=_DIG_ATTRS["curvature"])
        self.self.data_write("momlevel", description=_DIG_ATTRS["momlevel"])



class PInitInputData(clawpack.clawutil.data.ClawData):
    r"""
    D-Claw pressure initialization data object
    """

    def __init__(self, ndim):
        super(PInitInputData, self).__init__()

        # Set default values:
        self.add_attribute("init_ptype", 0)
        self.add_attribute("init_pmax_ratio", 1.0)
        self.add_attribute("init_ptf", 1.0)
        self.add_attribute("init_ptf2", 0.0)

    def write(self,out_file='setpinit.data',data_source='setrun.py'):
        self.open_data_file(out_file,data_source)   

        print("Creating data file setpinit.data")
        # open file and write a warning header:
        self.self.data_write("init_ptype", description=_DIG_ATTRS["init_ptype"])
        self.self.data_write("init_pmax_ratio", description=_DIG_ATTRS["init_pmax_ratio"])
        self.self.data_write("init_ptf", description=_DIG_ATTRS["init_ptf"])
        self.self.data_write("init_ptf2", description=_DIG_ATTRS["init_ptf2"])
