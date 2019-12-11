# nice: http://wiki.blender.org/index.php/Dev:2.5/Py/Scripts/Cookbook/Code_snippets/Interface
bl_info = { "name" : "Mount Electrodes", "category" : "Object" }

import os
import bpy
import bmesh
import mathutils
import sys
import os.path # see: http://stackoverflow.com/a/8384788
import math
import re      # regular expressions, to identify all the center-line positions *z

from mathutils import Vector
from math import sqrt

USE_DUMMY_ELECTRODE_POS=True

#
# introduce plugin to blender and setup scene properties
#
def register():
    bpy.utils.register_module( __name__ )

    # parse args for input path
    ARGV_SURFPATH_LBL      = "surfpath"
    ARGV_COORDINATE_LBL    = "coords"
    ARGV_SIZE_LBL          = "sizes"
    ARGV_MATRIX_LBL        = "transformmat"

    sinSurfPath    = ""
    sinSurfDir     = ""
    sinSurfFile    = ""
    coords         = ""
    matrixString   = "1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1"
    skinSurfPath   = ""
    sizes          = ""

    argv = sys.argv
    try:
        argv = argv[ argv.index( "--" ) + 1:]
    except ValueError:
        print( "No parameters provided, using default values!" )
        print( "Commandline argument options are: ")
        print( ARGV_SURFPATH_LBL+"\t...path to skin surface" )
        print( ARGV_COORDINATE_LBL+"\t\t...comma sepparated list of electrode coordinates (one 10/20 coordinate for each electrode)" )
        print( ARGV_SIZE_LBL+"\t\t...comma sepparated list of electrode sizes in the mm of the format NNxMM, eg. 50x75 for a 5cm x 10 cm electrode" )
        print( "Example: blender -- coords=Fz,Oz sizes=50x50x5,50x100x2 surfpath=/home/user/skin.stl" )
        argv = []

    for s in range( 0, len( argv) ):
        index = argv[s].find( "=" )
        if( index > -1 ):
            prefix=argv[s][:index]

            if( prefix == ARGV_SURFPATH_LBL ):
                skinSurfPath = argv[s][index+1:]
            elif( prefix == ARGV_COORDINATE_LBL ):
                coords = argv[s][index+1:]
            elif( prefix == ARGV_SIZE_LBL ):
                sizes = argv[s][index+1:]
            elif( prefix == ARGV_MATRIX_LBL ):
                matrixString = argv[s][index+1:]

    print( "+++++++ Parameters +++++++")
    print( "Path to skin-surface: " + skinSurfPath )
    print( "Coordinates of electrodes in 10/20 system: " + coords )
    print( "Size of electrode: " + sizes )
    print( "Matrix-String: " + matrixString )

    # process commandline arguments
    electrodeSizes = []
    electrodePositions = []
    if( len( sizes) > 0 ):
        electrodeSizes = [ pair.split('x') for pair in sizes.split(',') ]
    if( len( coords ) > 0 ):
        electrodePositions = coords.split(',')

    bpy.types.Scene.PredefinedElectrodeProperties = []
    for i in range (0,len(electrodePositions)):
        sizeCurrentElectrode = [ 50.0, 50.0, 5 ]
        if( i < len(electrodeSizes ) ):
                sizeCurrentElectrode = [ float(val) for val in electrodeSizes[i] ]
        bpy.types.Scene.PredefinedElectrodeProperties.append(                  \
            { 'size' : sizeCurrentElectrode, 'coord' : electrodePositions[i] } \
        )

    lenX = OBJECT_OT_MountElectrodes.DEFAULT_ELECTRODE_SIZE_X
    lenY = OBJECT_OT_MountElectrodes.DEFAULT_ELECTRODE_SIZE_Y
    lenZ = OBJECT_OT_MountElectrodes.DEFAULT_ELECTRODE_HEIGHT
    pos  = OBJECT_OT_MountElectrodes.DEFAULT_ELECTRODE_POSITION

    if len(bpy.types.Scene.PredefinedElectrodeProperties) > 0:
        lenX = bpy.types.Scene.PredefinedElectrodeProperties[0]['size'][0]
        lenY = bpy.types.Scene.PredefinedElectrodeProperties[0]['size'][1]
        lenY = bpy.types.Scene.PredefinedElectrodeProperties[0]['size'][2]
        pos   = bpy.types.Scene.PredefinedElectrodeProperties[0]['coord']

    # setup scene properties
    bpy.types.Scene.SkinSurfPath = bpy.props.StringProperty (
      name = "Path to surface file",
      default = skinSurfPath,
      description = "Specifiy file path to skin surface",
      subtype = 'FILE_PATH'
    )

    bpy.types.Scene.ElectrodePos = bpy.props.StringProperty (
      name = "Electrode position",
      default = pos,
      description = "Specify the position of Electrode in 10/20 coordinate system",
      subtype = "NONE"
    )

    bpy.types.Scene.MatrixString = bpy.props.StringProperty (
      name = "Transformation matrix",
      default = matrixString,
      description = "Specify the matrix used to transform the surface file into MNI space",
      subtype = "NONE"
    )

    bpy.types.Scene.ElectrodeX= bpy.props.FloatProperty (
      name = "Electrode length",
      default = lenX,
      description = "Specify the length (x) of Electrode in mm",
      min=1e-5,
      soft_min=1e-5,
      subtype="UNSIGNED"
    )

    bpy.types.Scene.ElectrodeY= bpy.props.FloatProperty (
      name = "Electrode width",
      default = lenY,
      description = "Specify the width (y) of Electrode in mm",
      min=1e-5,
      soft_min=1e-5,
      subtype="UNSIGNED"
    )

    bpy.types.Scene.ElectrodeZ= bpy.props.IntProperty (
      name = "Electrode height",
      default = lenZ,
      description = "Specify the height (z) of Electrode in mm",
      min=0,
      soft_min=0,
      subtype="UNSIGNED"
    )

    bpy.types.Scene.IDCurrentElectrode = bpy.props.IntProperty (
      name = "Electrode ID",
      default = 0,
      description = "Specify the electrode's ID (note: this property auto-increments)",
      min=0,
      soft_min=0,
      subtype="UNSIGNED"
    )

    bpy.types.Scene.ReloadEntity = bpy.props.EnumProperty (
        name        = "entity",
        description = "The enitity to reload",
        items       = [
            ("SKIN", "Skin", "", 0),
            ("SMOOTHSKIN", "Smooth Skin", "", 1),
            ("NILINE", "NI-Line", "", 2),
            ("LRLINE", "LR-Line", "", 3),
            ("HELPER", "Helper", "", 4)
        ]
    )

    bpy.types.Scene.AxisInverted = bpy.props.BoolProperty (
        name        = "Main axis inverted",
        default     = False,
        description = "Indiciate whether the main-axis (the y-axis) is aligned with or contrary to the subject's viewing axis"
    )
    bpy.types.Scene.UseSmootSkin = bpy.props.BoolProperty (
        name        = "Use smoot skin",
        default     = True,
        description = "Indiciate whether the smoothed version of the skin should be used for electrode modeling or the provided, original surface"
    )

    bpy.types.Scene.UseGelLayer = bpy.props.BoolProperty (
        name        = "Use gel layer",
        default     = True,
        description = "Create an additional gel layer between the electrode and the skin surface"
    )

    bpy.types.Scene.MyObjects = {}

#
# remove plugin from blenders list of active plugins
# remove scene properties
#
def unregister():
    bpy.utils.unregister_module( __name__ )
    bpy.context.scene.SkinSurfPath       = ""
    bpy.context.scene.ElectrodePos       = "Cz"
    bpy.context.scene.ElectrodeX         = 50.0
    bpy.context.scene.ElectrodeY         = 50.0
    bpy.context.scene.ElectrodeZ         = 5
    bpy.context.scene.IDCurrentElectrode = 0
    bpy.context.scene.ReloadEntity       = "SKIN"
    bpy.context.scene.AxisInverted       = False
    bpy.context.scene.UseSmootSkin       = True
    bpy.context.scene.UseGelLayer        = True
    bpy.context.scene.MyObjects.clear()

    del bpy.types.Scene.SkinSurfPath
    del bpy.types.Scene.ElectrodePos
    del bpy.types.Scene.ElectrodeX
    del bpy.types.Scene.ElectrodeY
    del bpy.types.Scene.ElectrodeZ
    del bpy.types.Scene.IDCurrentElectrode
    del bpy.types.Scene.ReloadEntity
    del bpy.types.Scene.AxisInverted
    del bpy.types.Scene.UseSmootSkin
    del bpy.types.Scene.UseGelLayer
    del bpy.types.Scene.MyObjects

class Utils( ):
    @staticmethod
    def parseMatrixString( _matS ):
        matrixElements = _matS.split(',')

        ret = mathutils.Matrix()

        index = 0
        for y in range(0,4):
            for x in range(0,4):
                ret[y][x] = float( matrixElements[index] )
                index += 1

        return ret

    @staticmethod
    def applyTransformation( _surface, _transformMat ):
        _surface.location = _transformMat * _surface.location

    @staticmethod
    def applyTransformationToVerts( _surface, _transformMat ):
        mesh = _surface.data
        for vertex in mesh.vertices:
            vertex.co = _transformMat * vertex.co

    @staticmethod
    def resetSelectedAndActiveObjects( ):
        if len( bpy.context.selected_objects ) > 0:
            bpy.ops.object.mode_set(mode='OBJECT')
            bpy.ops.object.select_all(action='DESELECT')
        bpy.context.scene.objects.active = None

    @staticmethod
    def selectAndActivateObject( _obj ):
        _obj.select = True
        bpy.context.scene.objects.active = _obj

    @staticmethod
    def showSelectedAndActive( ):
        print( "Active object" + repr( bpy.context.active_object ) )
        selected_objects = [ o for o in bpy.context.scene.objects if o.select ]
        print( "Selected objects: " + repr( selected_objects ) )
    
    #
    # Computes the center coordinate of a given object.
    # Assumptions:
    #   _electrode ... a bpy_types.Object (e.g. from bpy.types.Scene.MyObjects )
    @staticmethod
    def calcGeometricalCenter( _object ):
        obj_mesh = _object.to_mesh(bpy.context.scene, False, 'PREVIEW')
        center_x = [ v.co[0] for v in obj_mesh.vertices ]
        center_y = [ v.co[1] for v in obj_mesh.vertices ]
        center_z = [ v.co[2] for v in obj_mesh.vertices ]
        center = ( sum(center_x) / len(obj_mesh.vertices), \
                   sum(center_y) / len(obj_mesh.vertices), \
                   sum(center_z) / len(obj_mesh.vertices) )

        return center

    @staticmethod
    def calcMeshStatistics( _mesh ):
        # calc bounds
        minV = [ sys.maxsize, sys.maxsize, sys.maxsize ]
        maxV = [ 0, 0, 0 ]

        for point in range( 0, 8 ):
            for vectorComponent in range( 0, 3 ):
                minV[ vectorComponent ] =                               \
                    min(                                                \
                        _mesh.bound_box[ point ][ vectorComponent ],  \
                        minV[ vectorComponent ]                         \
                    )

                maxV[ vectorComponent ] =                               \
                    max(                                                \
                        _mesh.bound_box[ point ][ vectorComponent ],  \
                        maxV[ vectorComponent ]                         \
                    )

        bounds = { 'max' : Vector(maxV), 'min' : Vector(minV) }

        # calc center
        center = Utils.calcGeometricalCenter( _mesh )

        # return bounds, center and additional information
        return {
                    'bounds'      : bounds,
                    'center'      : Vector(center),
                    'origin'      : _mesh.location,
                    'dimensions'  : _mesh.dimensions
                }

    #
    # creates a one-vertex-object at the center of the provided mesh
    # useful for orienting other objects to the center of the provided mesh
    # @param _mesh  - the mesh the helper object should be created for
    #
    @staticmethod
    def createHelperObject( _mesh ):
        stats = Utils.calcMeshStatistics( _mesh )
        helperName = "helper_" + _mesh.name
        mesh   = bpy.data.meshes.new( helperName )
        helper = bpy.data.objects.new( helperName, mesh )
        helper.show_name = True

        # Link object to scene and save in custom array
        bpy.context.scene.objects.link( helper )

        # add vertices, create edges & faces and update
        mesh.from_pydata( [ stats['center'] ], [], [])
        mesh.update()
        Utils.resetSelectedAndActiveObjects()
        Utils.selectAndActivateObject( helper )
        bpy.ops.object.mode_set(mode='OBJECT')
        bpy.ops.object.origin_set( type='ORIGIN_GEOMETRY' )

        return helper
    #
    # creates new plane between the provided points with its center
    # at the center of the line connecting the points.
    # Furthermore the plane will aligned according to the thir point '_alignedTo'
    # @param _points    - an 2-element array with elements of type 'Vector'
    # @param _name      - a string naming the object to be created
    # @param _alignedTo - a 'Vector' denoting the point the plane should be tilted to
    #
    @staticmethod
    def createPlaneBetweenPoints( _points, _name, _alignedTo ):
        # create new object and insert the two points
        mesh = bpy.data.meshes.new( _name )
        obj  = bpy.data.objects.new( _name, mesh )
        obj.show_name = True

        # Link object to scene and save in custom array
        bpy.context.scene.objects.link( obj )
        bpy.types.Scene.MyObjects[ _name ] = obj

        # add vertices, create edge and update
        mesh.from_pydata( _points, [(0,1)], [])
        mesh.update()

        #extrude vertices in both sides
            # create extrusion vector, from center to twice the highest central point
        center = _points[0] - _points[1]
        center = _points[1] + 0.5 * center

        alignmentVector = _alignedTo - center

        Utils.resetSelectedAndActiveObjects()
        Utils.selectAndActivateObject( obj )
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.extrude_region_move( MESH_OT_extrude_region = { "mirror" : True },              \
                                          TRANSFORM_OT_translate = { "value"  : 2.0*alignmentVector } )
            # move the object half the way down to recenter it
        moveTo = -alignmentVector
        Utils.applyTransformation(              \
                obj,                            \
                Utils.parseMatrixString(        \
                "1,0,0,"  + repr( moveTo.x ) +  \
                ",0,1,0," + repr( moveTo.y ) +  \
                ",0,0,1," + repr( moveTo.z ) +  \
                ",0,0,0,1"                      \
                )                               \
        )
        bpy.ops.object.mode_set(mode='OBJECT')
        bpy.ops.object.origin_set(type='ORIGIN_GEOMETRY')

    #
    # assumptions on the state:
    # nothing in particular
    #
    @staticmethod
    def pojectFrameOntoSurface( _surf, _frame ):
        Utils.resetSelectedAndActiveObjects()
        Utils.selectAndActivateObject( _frame )

        bpy.ops.object.mode_set(mode='EDIT')
        planeMesh = bmesh.from_edit_mesh( _frame.data )

        # subdivide the edges so the line strip approximates
        # the surface in more detail
        for planeEdge in planeMesh.edges:
            planeEdge.select = True

        bpy.ops.mesh.subdivide( number_cuts=100 )
        bpy.ops.mesh.subdivide( number_cuts=100 )
       # bpy.ops.mesh.subdivide( number_cuts=2 )

        bpy.ops.object.mode_set(mode='OBJECT')

        # project the frame onto the skin-surface
        bpy.ops.object.modifier_add( type='SHRINKWRAP' )
        swm             = _frame.modifiers[0]
        swm.target      = _surf
        swm.wrap_method = "PROJECT"
        swm.use_negative_direction = True
        swm.use_positive_direction = True  # if set to False, the line will cross below plateaus
        swm.use_project_x = False
        swm.use_project_y = False
        swm.use_project_z = False

        bpy.ops.object.modifier_apply( apply_as='DATA', modifier=swm.name )

    @staticmethod
    def calcLineStats( _line ):
        Utils.resetSelectedAndActiveObjects()
        Utils.selectAndActivateObject( _line )

        bpy.ops.object.mode_set(mode='EDIT')
        lineMesh = bmesh.from_edit_mesh( _line.data )

        if not 'line_lengths' in bpy.types.Scene.MyObjects:
            bpy.types.Scene.MyObjects[ 'line_lengths' ] = {}

        bpy.types.Scene.MyObjects[ 'line_lengths' ][ _line.name ] = { \
            'totalLength' : 0.0,                                      \
            'edgesLength' : [None] * len( lineMesh.edges )            \
        }
          
        mat     = _line.matrix_world
        mat_inv = _line.matrix_world.inverted()

        # to measure the real leght in world coordinates
        Utils.applyTransformationToVerts( _line, mat )

        for ctr, edge in  enumerate( lineMesh.edges     ):
            # from blender-doc 'index':
            # This value is not necessarily valid, while editing the
            # mesh it can become dirty. It’s also possible to assign
            # any number to this attribute for a scripts internal logic.
            edge.index = ctr
            bpy.types.Scene.MyObjects[ 'line_lengths' ][ _line.name ][ 'edgesLength' ][ ctr ] =  \
            abs( edge.calc_length() )
            bpy.types.Scene.MyObjects[ 'line_lengths' ][ _line.name ][ 'totalLength' ]       +=  \
            bpy.types.Scene.MyObjects[ 'line_lengths' ][ _line.name ][ 'edgesLength' ][ ctr ]

        Utils.applyTransformationToVerts( _line, mat_inv )

        print( "Total length of line '" + _line.name + "' " + repr(bpy.types.Scene.MyObjects[ 'line_lengths' ][ _line.name ][ 'totalLength' ]) + " (object space!)" )


    #
    # assumptions on the state of '_mesh'
    # - is a BMesh (mesh in edit mode)
    #
    @staticmethod
    def deselectAllVertices( _mesh ):
        assert type( _mesh ) == bmesh.types.BMesh

        for vertex in _mesh.verts:
            vertex.select = False

    @staticmethod
    def deselectAllEdges( _mesh ):
        assert type( _mesh ) == bmesh.types.BMesh

        for edge in _mesh.edges:
            edge.select = False

    @staticmethod
    def deselectAllFaces( _mesh ):
        assert type( _mesh ) == bmesh.types.BMesh

        for face in _mesh.faces:
            face.select = False

    #
    # assumption on the state of '_mesh'
    # - is an active BMesh ( mesh in edit mode )
    # Note: any reference to the returned 'lowestVertex' will become
    #       invalid, if you leave EDIT-mode of '_mesh'
    #
    @staticmethod
    def findLowestVertex( _mesh, _spatialConstraint = lambda _vertex : True ):
        assert type( _mesh ) == bmesh.types.BMesh

        if hasattr(_mesh.verts, "ensure_lookup_table"):
            _mesh.verts.ensure_lookup_table()

        lowestVertex = _mesh.verts[0]
        for vertex in _mesh.verts:
            vertex.select = False
            if vertex.co.z < lowestVertex.co.z and  \
               _spatialConstraint( vertex.co ) :
                lowestVertex = vertex

        return lowestVertex


    #
    # assumptions:
    # @param _start     is on of the two endings of '_strip' and the
    #                   reference point for '_fraction'
    # @param _strip     is of type 'bmesh.typesBMesh' (bmesh in edit mode)
    # @param _lineStats use precalculated line-statistics for length measurment
    #                   (more precise)
    # @param _direction denotes which axis to use to compare consecutive vertices
    # @param _forward   which direction to go: declining or increasing along '_direction'
    #
    @staticmethod
    def vertexAtLinestripFraction( _strip, _start, _fraction, _lineStats=None, _direction='z', _forward=True  ):
        cmpr = None

        if _forward:
            cmpr = lambda fl1, fl2 : fl1 > fl2
        else:
            cmpr = lambda fl1, fl2 : fl1 < fl2

        # dtermine the vector component to compare
        dim = 0
        if _direction.lower() == 'x' :
            dim = 0
        elif _direction.lower() == 'y' :
            dim = 1
        elif _direction.lower() == 'z' :
            dim = 2
        else:
            print( "Invalid _direction value provided! ('x', 'y', 'z' are accepted) ")
            return None

        # determine which vertex to follow, which edge to go (startEdge)
        # (if the start-vertex is in the middle of the line
        #  and thus has two neighbours)
        currentVert   = _start
        lastVisit     = _start
        nextVert      = None

        startEdges = _start.link_edges
        startEdge  = None

        for edge in startEdges:
            nextVert = edge.other_vert( _start )
            if cmpr( nextVert.co[ dim ], _start.co[ dim ] ):
                currentVert = nextVert
                startEdge   = edge
                break

        if startEdge is None:
            print( "Trasversal direction is ambiguous! Cannot traverse." )
            return None

        # traverse the whole line from the beginning if no
        # line-statistics provided until we reach the vertex
        # at '_fraction'
        if _lineStats is None:
            totalNumVerts = len( _strip.verts ) - 1
            for v in range( 1, int( math.ceil( totalNumVerts * _fraction ) ) ):
                edges = currentVert.link_edges
                for edge in edges:
                    nextVert = edge.other_vert( currentVert )
                    if not nextVert == lastVisit:
                        lastVisit   = currentVert
                        currentVert = nextVert
                        break
        # use the line statistics to calculate the vertex at '_fraction'
        else:
            maxLength     = _lineStats['totalLength'] * _fraction
            currentLenght = _lineStats['edgesLength'][ startEdge.index ]
            while currentLenght < maxLength:
                edges = currentVert.link_edges
                for edge in edges:
                    nextVert = edge.other_vert( currentVert )
                    if not nextVert == lastVisit:
                        lastVisit   = currentVert
                        currentVert = nextVert
                        currentLenght += _lineStats['edgesLength'][ edge.index ]
                        break


        return currentVert

    @staticmethod
    def euklidianDist( _pt1, _pt2 ):
        xDiff = _pt1.x - _pt2.x
        yDiff = _pt1.y - _pt2.y
        zDiff = _pt1.z - _pt2.z

        return sqrt( xDiff*xDiff + yDiff*yDiff + zDiff*zDiff )
    #
    # assumptions on the state:
    # @param _mesh          is a BMesh in edit mode
    # @param _referencePt   is a 'Vector' representing the query-point
    #
    @staticmethod
    def findClosestMeshPoint( _object, _referencePt ):
        Utils.resetSelectedAndActiveObjects( )
        Utils.selectAndActivateObject( _object )

        bpy.ops.object.mode_set(mode='EDIT')
        mesh = bmesh.from_edit_mesh( _object.data )

        refPtLocal = _object.matrix_world.inverted() * _referencePt

        nearest = None
        minDist = sys.maxsize
        dist    = 0
        for vertex in mesh.verts:
            dist = Utils.euklidianDist( vertex.co, refPtLocal )
            if dist < minDist:
                minDist = dist
                nearest = vertex

        return {'point' : nearest, 'distance' : minDist, 'mesh' : mesh }


    @staticmethod
    def removeVerticesAccordingToDirection( _meshOfObject, _angle, _direction ):
        bmeshOfObject = bmesh.from_edit_mesh( _meshOfObject.data )
        obj = { 'mesh' : _meshOfObject, 'bmesh' : bmeshOfObject, 'components' : bmeshOfObject.verts }
        Utils.removeMeshPartsAccordingToDirection( obj, _angle, _direction, 1 )

    @staticmethod
    def removeFacesAccordingToDirection( _meshOfObject, _angle, _direction, ):
        bmeshOfObject = bmesh.from_edit_mesh( _meshOfObject.data )
        obj = { 'mesh' : _meshOfObject, 'bmesh' : bmeshOfObject, 'components' : bmeshOfObject.faces }
        Utils.removeMeshPartsAccordingToDirection( obj, _angle, _direction, 5 )

    @staticmethod
    def removeMeshPartsAccordingToDirection( _obj, _angle, _direction, _meshComponent ):
            # vertex normals are automatically calculated by Blender each time one
            # toggles between object-mode and edit mode.
            # How vertex normals are calulcated:
            #  > If the vertex pertains to one or more faces, its normal is
            #    the interpolated value of all the normals of these faces.
            #  > If the vertex does not pertain to any face, its normal is
            #    aligned with the line from the object’s center passing
            #    through this vertex.
        bpy.ops.mesh.normals_make_consistent()

        componentsToRemove= [ ]
        for component in _obj['components']:
            nLocal   = component.normal.to_4d()
            nLocal.w = 0
            nWorld   = ( _obj['mesh'].matrix_world * nLocal ).to_3d()
            nWorld.normalize()

            dot = nWorld.dot( _direction )

                # remove all faces with an angle smaller than ~108 degree to 'direction'
                # 'direction' points to the skin surface, the normals of the face we want to
                # keep point outwards the skin surface
                # cos(x) < 0 = angles between 90degreed and 270 degrees
            if dot > _angle:
                componentsToRemove.append( component )

                # the integer for 'context' in 'delete' has to be one of the values of the enum defined in
                # 'bmesh_operator_api.h' :
                # DEL_VERTS = 1,DEL_EDGES,DEL_ONLYFACES,DEL_EDGESFACES,DEL_FACES,DEL_ALL,DEL_ONLYTAGGED
        bmesh.ops.delete( _obj['bmesh'], geom=componentsToRemove, context=_meshComponent)
        bmesh.update_edit_mesh( _obj['mesh'].data, True) # make changes to the mesh persitent, recalculate normals

        Utils.deselectAllFaces( _obj['bmesh'] )

    @staticmethod
    def populateInputs( _context, _electrodeID ):
        if( len( _context.scene.PredefinedElectrodeProperties ) > 0 ):
            _context.scene.ElectrodePos = _context.scene.PredefinedElectrodeProperties[ _electrodeID ]['coord']
            _context.scene.ElectrodeX   = _context.scene.PredefinedElectrodeProperties[ _electrodeID ]['size'][0]
            _context.scene.ElectrodeY   = _context.scene.PredefinedElectrodeProperties[ _electrodeID ]['size'][1]
            _context.scene.ElectrodeZ   = _context.scene.PredefinedElectrodeProperties[ _electrodeID ]['size'][2]
        else:
            _context.scene.ElectrodePos = OBJECT_OT_MountElectrodes.DEFAULT_ELECTRODE_POSITION
            _context.scene.ElectrodeX   = OBJECT_OT_MountElectrodes.DEFAULT_ELECTRODE_SIZE_X
            _context.scene.ElectrodeY   = OBJECT_OT_MountElectrodes.DEFAULT_ELECTRODE_SIZE_Y
            _context.scene.ElectrodeY   = OBJECT_OT_MountElectrodes.DEFAULT_ELECTRODE_HEIGHT

    @staticmethod
    def resetAndEnableInputs( _context ):
        _context.scene.ElectrodePos =  OBJECT_OT_MountElectrodes.DEFAULT_ELECTRODE_POSITION
        _context.scene.ElectrodeX   =  OBJECT_OT_MountElectrodes.DEFAULT_ELECTRODE_SIZE_X
        _context.scene.ElectrodeY   =  OBJECT_OT_MountElectrodes.DEFAULT_ELECTRODE_SIZE_Y
        _context.scene.ElectrodeZ   =  OBJECT_OT_MountElectrodes.DEFAULT_ELECTRODE_HEIGHT

    #
    # Uses the helper grid to get the world coordinates of the 10/20 coordinate
    # Assumptions:
    #   _1020_coord ... string ... Identifier of the 1020 corrdinate (e.g. C3) 
    #
    @staticmethod
    def obtainElectrodeCoord_from1020( _1020coord ):
        coord = bpy.context.scene.MyObjects['1020grid'].getWorldCoord( _1020coord )

        if coord:
            return coord
        else:
            return Vector( ( 23.0, 80.0, 220.0 ) )
    #
    # Uses the electrode-center to calculate the world corrdinates of the electrode
    # Note: other than the 'obtainElectrodeCoord_from1020' method this coordinate is
    #       not immediately situated on the skin surface.
    #
    @staticmethod
    def obtainElectrodeCoord_center( _electrode ):
        mat     = _electrode.matrix_world
        mat_inv = _electrode.matrix_world.inverted()

        Utils.applyTransformationToVerts( _electrode, mat )

        center = Utils.calcGeometricalCenter( _electrode ) 
        
        Utils.applyTransformationToVerts( _electrode, mat_inv )
        
        # to measure the real leght in world coordinates
        return Vector( center )


#
# popup which will be displayed, when an error occurs
#
class OBJECT_OT_PopUp( bpy.types.Operator ):
    bl_idname   = "error.message"
    bl_label    = "Message"

    mType       = bpy.props.StringProperty()
    mMessage    = bpy.props.StringProperty()

    def execute(self, context):
        self.report({'INFO'}, self.mMessage)
        print(self.mMessage)
        return {'FINISHED'}

    def invoke(self, context, event):
        windowmgr = context.window_manager
        return windowmgr.invoke_popup(self, width=768, height=200)

    def draw(self, context):
        self.layout.label("MountElectrodes-output: ")
        row = self.layout.split(0.15)
        row.prop(self, "mType")
        row.prop(self, "mMessage")

#
# helper operator, that loads the selected skin surface
#
class OBJECT_OT_LoadSkinSurf( bpy.types.Operator ):
    bl_idname       = "object.load_skin_surface"
    bl_label        = "Load skin surface"
    bl_options      = { "REGISTER", "UNDO" }
    bl_description  = "Loads the skin surface, the electordes will be mounted onto"

    def execute( _self, _context ):
        # parse the path from the scene string property
        skinSurfFile = os.path.basename( _context.scene.SkinSurfPath )
        skinSurfDir  = os.path.dirname( _context.scene.SkinSurfPath )

        # import
        imported = bpy.ops.import_mesh.stl  \
        (                                   \
            filepath=skinSurfFile,          \
            filter_glob="*.stl",            \
            files=[{"name":skinSurfFile}],  \
            directory=skinSurfDir           \
        )

        # skin surface is automatically selected after import, save the reference to that object
        selectedObjects = [ o for o in bpy.context.scene.objects if o.select ]
        bpy.types.Scene.MyObjects["skin"] = selectedObjects[0]
            # we have write-access to scene-variables through bpy.types.Scene

        # transform surface into MNI space
        #~ bpy.ops.object.origin_set(type='GEOMETRY_ORIGIN')
        Utils.applyTransformation(                                      \
                _context.scene.MyObjects["skin"],                       \
                Utils.parseMatrixString( _context.scene.MatrixString )  \
        )
            # we have read access to scene-variables trough the current context
            
        # create a smoothed version of the skin; we will use it later to obtain smooth normals from
        bpy.ops.object.duplicate()
        smoothSkinSurface = bpy.context.active_object
            # we use the laplacian smoothing algorithm as it is more feature preserving
        bpy.ops.object.modifier_add(type='LAPLACIANSMOOTH')
        laplacianSmooth               = smoothSkinSurface.modifiers[0]
        laplacianSmooth.lambda_factor = 2
        laplacianSmooth.iterations    = 30
        
           # in order to apply a modifier two preconditions have to satisfied:
           # a) the object whose modifiers should be applied needs to be the active object
           # b) the name of the modifier to be applied has to be provided (e.g. csg.name)
        bpy.ops.object.modifier_apply( apply_as='DATA', modifier=laplacianSmooth.name )
            
        bpy.types.Scene.MyObjects[ "smoothskin" ] = smoothSkinSurface
        
        # turn viewport to surface
        bpy.ops.view3d.view_selected()
        bpy.ops.view3d.view_persportho()

        # precalculate mesh statistics
        bpy.types.Scene.MyObjects["skin_stats"] = Utils.calcMeshStatistics( bpy.types.Scene.MyObjects["skin"] )
        bpy.types.Scene.MyObjects["helper"]     = Utils.createHelperObject( bpy.types.Scene.MyObjects["skin"] )

        #success, always...
        return {"FINISHED"}

#
# helper operator to export the created elextrode surfaces as STL files
#
class OBJECT_OT_ExportElectrodes( bpy.types.Operator ):
    """Export all electrode parts of scene as surfaces"""

    bl_idname       = "object.export_electrodes"
    bl_label        = "Export Electrodes"
    bl_options      = { "REGISTER", "UNDO" }
    bl_description  = "Exports all electrode-part of the scene"

    mPath = bpy.props.StringProperty( name="Output-Path",                               \
                                      description="Define the location where to store the electrodes", \
                                      default="/tmp/" )

    def invoke( _self, _context, _event):
        windowmgr = _context.window_manager
        return windowmgr.invoke_props_dialog( _self )

    def execute( _self, _context ):
        scene   = _context.scene

        electrodeParts = [ part for part in scene.objects if part.name.startswith("electrode") ]

        for p in range(0, len( electrodeParts ) ):
            for obj in bpy.context.selectable_objects:
                obj.select = False
            bpy.ops.object.select_pattern( pattern=electrodeParts[ p ].name, case_sensitive=False )
            bpy.context.scene.objects.active = electrodeParts[ p ]
            print( "Exporting: " + os.path.join(_self.mPath, electrodeParts[ p ].name + ".stl" )  )
            bpy.ops.export_mesh.stl                                                     \
            (                                                                           \
                filepath=os.path.join(_self.mPath, electrodeParts[ p ].name + ".stl" ), \
                ascii=True,                                                             \
                use_mesh_modifiers=True,                                                \
                check_existing=False,                                                   \
                use_selection=True                                                      \
            )

        return {"FINISHED"}

class OBJECT_OT_ReloadSkinSurface( bpy.types.Operator ):
    """Reset internal data-structure without actually lading the surface from a STL file"""

    bl_idname       = "object.reload_skinsurf"
    bl_label        = "Reload Skin Surface"
    bl_options      = { "REGISTER", "UNDO" }
    bl_description  = "Reloads the skin surface"

    def execute( _self, _context ):
        # assume that the skin surface is the only currently selected object
        selectedObjects = [ o for o in bpy.context.scene.objects if o.select ]

        if( _context.scene.ReloadEntity == "SKIN" ):
            bpy.types.Scene.MyObjects["skin"] = selectedObjects[0]
            # precalculate mesh statistics
            bpy.types.Scene.MyObjects["skin_stats"] = Utils.calcMeshStatistics( bpy.types.Scene.MyObjects["skin"] )
        elif( _context.scene.ReloadEntity == "NILINE" ):
            bpy.types.Scene.MyObjects[ OBJECT_OT_CreateNasionInionPlane.PLANE_NAME ] = selectedObjects[0]
        elif( _context.scene.ReloadEntity == "LRLINE"):
            bpy.types.Scene.MyObjects[ OBJECT_OT_CreateLRPlane.PLANE_NAME ] = selectedObjects[0]
        elif( _context.scene.ReloadEntity == "HELPER"):
            bpy.types.Scene.MyObjects[ "helper" ] = selectedObjects[0]       
        elif( _context.scene.ReloadEntity == "SMOOTHSKIN"):
            bpy.types.Scene.MyObjects[ "smoothskin" ] = selectedObjects[0]

        return {"FINISHED"}

#
# layout code, for side-panel with buttons and text-input
#
class OBJECT_PT_MountElectrodes( bpy.types.Panel ):
    bl_label          = "Mount-Electrodes"
    bl_category       = "Mount Electrodes"
    bl_space_type     = "VIEW_3D"
    bl_region_type    = "TOOLS"

    def draw( _self, _context):
        _self.layout.prop( _context.scene, "ReloadEntity")
        _self.layout.operator( OBJECT_OT_ReloadSkinSurface.bl_idname, text='Reload')
        _self.layout.prop( _context.scene, "SkinSurfPath" )
        _self.layout.operator(OBJECT_OT_LoadSkinSurf.bl_idname, text='Load Surface')
        _self.layout.prop( _context.scene, "AxisInverted" )
        _self.layout.prop( _context.scene, "UseSmootSkin" )
        _self.layout.operator(OBJECT_OT_CreateNasionInionPlane.bl_idname, text='Set NI-Plane')
        _self.layout.operator(OBJECT_OT_CreateNasionInionHelperLine.bl_idname, text='Create NI-Helper-Line')
        _self.layout.operator(OBJECT_OT_SetNasion.bl_idname, text='Set Nasion')
        _self.layout.operator(OBJECT_OT_SetInion.bl_idname, text='Set Inion')
        _self.layout.operator(OBJECT_OT_CreateLRPlane.bl_idname, text='Set LR-Plane')
        _self.layout.operator(OBJECT_OT_CreateLRHelperLine.bl_idname, text='Create LR-Helper-Line')
        _self.layout.operator(OBJECT_OT_SetLPoint.bl_idname, text='Set L-Point')
        _self.layout.operator(OBJECT_OT_SetRPoint.bl_idname, text='Set R-Point')
        _self.layout.prop( _context.scene, "ElectrodePos" )
        _self.layout.prop( _context.scene, "ElectrodeX" )
        _self.layout.prop( _context.scene, "ElectrodeY" )
        _self.layout.prop( _context.scene, "ElectrodeZ" )
        _self.layout.prop( _context.scene, "IDCurrentElectrode" )
        _self.layout.operator(OBJECT_OT_ElectrodePrototype.bl_idname, text='Set Prototype')
        _self.layout.prop( _context.scene, "UseGelLayer" )
        _self.layout.operator(OBJECT_OT_MountElectrodes.bl_idname, text='Set electrode')

class OBJECT_OT_CreateNasionInionPlane( bpy.types.Operator ):
    """create the helper plane to define the nasion an inion"""

    bl_idname       = "object.create_niplane"
    bl_label        = "Create NI-plane"
    bl_options      = { "REGISTER", "UNDO" }
    bl_description  = "Erstellt die Hilfsebene zwischen Nasion und Inion"

    PLANE_NAME = "niLine"

    def execute( _self, _context ):
        modelStats =  _context.scene.MyObjects[ "skin_stats" ]

        leftP  = Vector( ( modelStats['center'].x, modelStats['bounds']['min'].y, modelStats['center'].z) )
        rightP = Vector( ( modelStats['center'].x, modelStats['bounds']['max'].y, modelStats['center'].z) )
        Utils.createPlaneBetweenPoints(    \
            [                              \
                leftP,                     \
                rightP                     \
            ],                             \
            _self.PLANE_NAME,              \
            Vector( (modelStats['center'].x, modelStats['center'].y, modelStats['bounds']['max'].z) )
        )

        _context.scene.MyObjects[ _self.PLANE_NAME ].lock_location = [ False, True, True ]

        return {"FINISHED"}

class OBJECT_OT_CreateLRPlane( bpy.types.Operator ):
    """create the helper plane between the ears"""

    bl_idname       = "object.create_lrplane"
    bl_label        = "Create LR-plane"
    bl_options      = { "REGISTER", "UNDO" }
    bl_description  = "Erstellt die Hilfsebene zwischen den Ohren"

    PLANE_NAME = "lrLine"

    def execute( _self, _context ):
        modelStats = _context.scene.MyObjects[ "skin_stats" ]
        niLine     = _context.scene.MyObjects[ OBJECT_OT_CreateNasionInionPlane.PLANE_NAME ]

        leftP  = Vector( ( modelStats['bounds']['min'].x, modelStats['center'].y, modelStats['center'].z) )
        rightP = Vector( ( modelStats['bounds']['max'].x, modelStats['center'].y, modelStats['center'].z) )

        # determine center of NI-line; will be used to align the lrPlane
        Utils.calcLineStats( bpy.context.scene.MyObjects[ OBJECT_OT_CreateNasionInionPlane.PLANE_NAME ] )

        Utils.resetSelectedAndActiveObjects()
        Utils.selectAndActivateObject( niLine )

        bpy.ops.object.mode_set(mode='EDIT')
        niLineMesh     = bmesh.from_edit_mesh( niLine.data )
        totalNumVerts  = len( niLineMesh.verts )
        lowestVertex   = Utils.findLowestVertex( niLineMesh )
        centerOfNILine = Utils.vertexAtLinestripFraction( niLineMesh, lowestVertex, 0.5, bpy.types.Scene.MyObjects[ 'line_lengths' ][ niLine.name ] )
        planeTarget    = niLine.matrix_world * centerOfNILine.co

        Utils.createPlaneBetweenPoints(    \
            [                              \
                leftP,                     \
                rightP                     \
            ],                             \
            _self.PLANE_NAME,              \
            planeTarget
        )

        _context.scene.MyObjects[ _self.PLANE_NAME ].lock_location = [ True, False, True ]

        return {"FINISHED"}

class SetMarker( bpy.types.Operator ):
    bl_idname  = "object.setmarker_super"
    bl_label   = "setMarker"
    bl_options = { "REGISTER", "UNDO" }

    mLine    = None
    fCompare = None

    def __init__( _self, _line, _compare ):
        _self.mLine     = _line
        _self.fCompare  = _compare
    #
    # assumptions about the state:
    # - '_line' is an object in edit mode
    # - '_lineMesh' is a BMesh from '_line' (mesh in edit mode)
    # - '_markerData' is a dictionary, containing
    #   * 'marker'  : the marker as a BMVert from '_lineMesh'
    #   * 'compare' : and additional lambda function which compares
    #                 the currently checked vertex with the model-center
    #                 ( for easier distinction between front/back, left/right )
    #
    @staticmethod
    def removeBelowSpecialMarker( _context, _line, _lineMesh, _markerData ):
        modelStats  = _context.scene.MyObjects[ "skin_stats" ]
        localCenter = _line.matrix_world.inverted() * modelStats['center']

        # check if the marker itself is on the correct side of the object
        if not _markerData[ 'compare']( _markerData[ 'marker' ].co, localCenter ):
            return

        print( _markerData[ 'compare'] )

        toBeRemoved = []
        for lineVertex in _lineMesh.verts:
            if lineVertex.co.z < _markerData['marker'].co.z and                   \
               _markerData['compare'](lineVertex.co, localCenter ) :
               toBeRemoved.append( lineVertex )
                                                           # 1 = vertices
        bmesh.ops.delete( _lineMesh, geom=toBeRemoved, context=1)
        bmesh.update_edit_mesh( _line.data )

    def execute( _self, _context ):
        modelStats = _context.scene.MyObjects[ "skin_stats" ]

        Utils.resetSelectedAndActiveObjects()
        Utils.selectAndActivateObject( _self.mLine )
        bpy.ops.object.mode_set(mode='EDIT')

        lineMesh = bmesh.from_edit_mesh( _self.mLine.data )

        # finde the selected vertex
        marker = None
        for lineVertex in lineMesh.verts:
            if lineVertex.select:
                marker = lineVertex
                break

        SetMarker.removeBelowSpecialMarker(
            _context,
            _self.mLine,
            lineMesh,
            { 'marker'  : marker,
              'compare' : _self.fCompare }
        )

        return {"FINISHED"}

#
# assumptions about the state:
# - 'skin_stats' and the nasion/inion-line exist
# - only one vertex of the nasion/inion-line is selected
#
class OBJECT_OT_SetNasion( SetMarker ):
    """set the nasion, must be slected already"""

    bl_idname       = "object.set_nasion"
    bl_label        = "Set nasion"
    bl_options      = { "REGISTER", "UNDO" }
    bl_description  = "Passt die Hilfsgerade auf das Nasion an"

    # We need to declare fCompare outside the constructor because it is
    # sometimes used as static member of 'OBJECT_OT_SetNasion'.
    # Also we cannot branch according to bpy.context.scene.AxisInverted here
    # because this porerty is not necessarily already set upon parsing
    # of OBJECT_OT_SetNasion (not construction, we are outside the
    # constructor here!!). So we put the decision for the branch within
    # the lambda. As lambda are not allowed to be multi-line commands
    # and there is no shorthand-if-then-else in python we need to branch
    # by clever logical operation:
    # The expression evaluates TRUE, if
    #  a) the nose is in the negative direction to the center
    #           AND
    #     the axis is inverted (person looking negatively down the y axis to the origin)
    #  b) the nose is in the positive direction to the center
    #           AND
    #     the axis is not inverted (person looking positively along with the y-axis towards the origin)
    fCompare = lambda vertex, center : (vertex.y < center.y) == bpy.context.scene.AxisInverted

    def __init__( _self ):
        super( OBJECT_OT_SetNasion, _self).__init__(
            bpy.context.scene.MyObjects[ OBJECT_OT_CreateNasionInionPlane.PLANE_NAME ],
            OBJECT_OT_SetNasion.fCompare
        )

    def execute( _self, _context ):
        return super( OBJECT_OT_SetNasion, _self).execute( _context )

# assumptions about the state:
# - 'skin_stats' and the nasion/inion-line exist
# - only one vertex of the nasion/inion-line is selected
#
class OBJECT_OT_SetInion( SetMarker ):
    """set the inion, must be slected already"""

    bl_idname       = "object.set_inion"
    bl_label        = "Set inion"
    bl_options      = { "REGISTER", "UNDO" }
    bl_description  = "Passt die Hilfsgerade auf das Inion an"

    fCompare = lambda vertex, center : (vertex.y > center.y) ==  bpy.context.scene.AxisInverted

    def __init__( _self ):
        super( OBJECT_OT_SetInion, _self).__init__(
            bpy.context.scene.MyObjects[ OBJECT_OT_CreateNasionInionPlane.PLANE_NAME ],
            OBJECT_OT_SetInion.fCompare
        )

    def execute( _self, _context ):
        return super( OBJECT_OT_SetInion, _self).execute( _context )

class OBJECT_OT_SetLPoint( SetMarker ):
    """set the l point, must be slected already"""

    bl_idname       = "object.set_lpoint"
    bl_label        = "Set L-point"
    bl_options      = { "REGISTER", "UNDO" }
    bl_description  = "Passt die Hilfsgerade auf das den L-Punkt an"


    def __init__( _self ):
        if bpy.context.scene.AxisInverted:
            fCompare = lambda vertex, center : vertex.x < center.x
        else:
            fCompare = lambda vertex, center : vertex.x > center.x

        super( OBJECT_OT_SetLPoint, _self).__init__(
            bpy.context.scene.MyObjects[ OBJECT_OT_CreateLRPlane.PLANE_NAME ],
            fCompare
        )

    def execute( _self, _context ):
        return super( OBJECT_OT_SetLPoint, _self).execute( _context )

class OBJECT_OT_SetRPoint( SetMarker ):
    """set the r point, must be slected already"""

    bl_idname       = "object.set_rpoint"
    bl_label        = "Set R-point"
    bl_options      = { "REGISTER", "UNDO" }
    bl_description  = "Passt die Hilfsgerade auf das den R-Punkt an"

    def __init__( _self ):
        if bpy.context.scene.AxisInverted:
            fCompare = lambda vertex, center : vertex.x > center.x
        else:
            fCompare = lambda vertex, center : vertex.x < center.x

        super( OBJECT_OT_SetRPoint, _self).__init__(
            bpy.context.scene.MyObjects[ OBJECT_OT_CreateLRPlane.PLANE_NAME ],
            fCompare
        )

    def execute( _self, _context ):
        return super( OBJECT_OT_SetRPoint, _self).execute( _context )

class CreateHelperLine( bpy.types.Operator ):
    bl_idname  = "object.helper_lone_super"
    bl_label   = "createHelperLine"
    bl_options = { "REGISTER", "UNDO" }

    mPlane = None

    def __init__( _self, _plane ):
        _self.mPlane = _plane

    def execute( _self, _context ):
        modelStats  = _context.scene.MyObjects[ "skin_stats" ]

        bpy.ops.object.mode_set(mode='OBJECT')

        Utils.resetSelectedAndActiveObjects()
        Utils.selectAndActivateObject( _self.mPlane )

        bpy.ops.object.mode_set(mode='EDIT')
        planeMesh = bmesh.from_edit_mesh( _self.mPlane.data )

        # clear any previous selection
        Utils.deselectAllVertices( planeMesh )
        Utils.deselectAllEdges( planeMesh )

        # select lowest horizontal edge
        toBeRemoved = []
        for planeEdge in planeMesh.edges:
            if planeEdge.verts[0].co.z == planeEdge.verts[1].co.z:
                globalCoord = _self.mPlane.matrix_world * planeEdge.verts[0].co
                if globalCoord.z < modelStats['center'].z:
                    planeEdge.select = True
                    toBeRemoved.append( planeEdge )
                    break

        # and remove it to create a frame around the head   # 4 = edgesfaces - don't know why not 2 = edges...
        bmesh.ops.delete(planeMesh, geom=toBeRemoved, context=4)
        bmesh.update_edit_mesh( _self.mPlane.data )

        # project the frame onto the head and subdivide it
        Utils.pojectFrameOntoSurface(            \
            _context.scene.MyObjects[ "skin" ],  \
            _self.mPlane                         \
         )

        return {"FINISHED"}

class OBJECT_OT_CreateLRHelperLine( CreateHelperLine ):
    """Create the helper line to specify the left and right point"""

    bl_idname       = "object.create_lrline"
    bl_label        = "Create LR-line"
    bl_options      = { "REGISTER", "UNDO" }
    bl_description  = "Erstellt Hilfslinie zur Auswahl des linken und rechten Punktes"

    def __init__( _self ):
        super( OBJECT_OT_CreateLRHelperLine, _self).__init__(               \
            bpy.context.scene.MyObjects[ OBJECT_OT_CreateLRPlane.PLANE_NAME ]  \
        )

    #
    # Assumptions on the state:
    # - helper plane exists
    # - head surface exists
    #
    def execute( _self, _context ):
        res = super( OBJECT_OT_CreateLRHelperLine, _self ).execute( _context )

        # enter edit mode for L/R point selection
        bpy.ops.object.mode_set(mode='EDIT')

        return res


class OBJECT_OT_CreateNasionInionHelperLine( CreateHelperLine ):
    """Create the helper line to specify the nasion and the inion later"""

    bl_idname       = "object.create_niline"
    bl_label        = "Create NI-line"
    bl_options      = { "REGISTER", "UNDO" }
    bl_description  = "Erstellt Hilfslinie zur Auswahl von Nasion & Inion"

    xyPlane = None

    def __init__( _self ):
        _self.xyPlane = bpy.context.scene.MyObjects[ OBJECT_OT_CreateNasionInionPlane.PLANE_NAME ]

        super( OBJECT_OT_CreateNasionInionHelperLine, _self).__init__(  \
            _self.xyPlane                                               \
        )

    #
    # Assumptions on the state:
    # - helper plane exists
    # - head surface exists
    #
    def execute( _self, _context ):
        super( OBJECT_OT_CreateNasionInionHelperLine, _self ).execute( _context )

        # try to find the nasion on the frame
            # find the tip of the nose
        bpy.ops.object.mode_set(mode='EDIT')
        lineMesh = bmesh.from_edit_mesh( _self.xyPlane.data )

        Utils.deselectAllVertices( lineMesh )
        Utils.deselectAllEdges( lineMesh )

            # need to call this before direct indexing is possible
        if hasattr(lineMesh.verts, "ensure_lookup_table"):
            lineMesh.verts.ensure_lookup_table()

        nosetip = lineMesh.verts[0]
        for lineVertex in lineMesh.verts:
            lineVertex.select = False          # take axis inversion into consideration! (see comment in OBJECT_OT_SetNasion for furhter details!)
            if (lineVertex.co.y < nosetip.co.y) == bpy.context.scene.AxisInverted:
                nosetip = lineVertex

            # the nasion is at the inflection-point of the helper line
        nasionNotFound = True
        abort          = False
        currentVert    = nosetip
        numInvalidPts  = 0
        while nasionNotFound and not numInvalidPts == 2:
            edges = currentVert.link_edges
            # climb upwards
            numInvalidPts = 0
            for edge in edges:
                nextVert = edge.other_vert( currentVert )
                if nextVert.co.z > currentVert.co.z:    # take axis inversion into consideration! (see comment in OBJECT_OT_SetNasion for furhter details!)
                    if (nextVert.co.y < currentVert.co.y) == bpy.context.scene.AxisInverted:
                        nasionNotFound = False
                        nasion = currentVert
                        print( "Nasion found at" )
                        print( nasion.co )
                    else:
                        currentVert = nextVert
                    break
                else:
                    numInvalidPts += 1 # a point is invalid when we climb down again

        if numInvalidPts == 2:
            bpy.ops.error.message('INVOKE_DEFAULT',
                                  mType = "Error",
                                  mMessage = "Nasion could not be found! Please specify manually.")
            print( "Could not find nasion!" )
            # deselect all vertices to allow the user to manually select one
            Utils.deselectAllVertices( lineMesh )
            return {"FINISHED"}

            # remove all vertices below the nasion on the front side of the head
        SetMarker.removeBelowSpecialMarker(                 \
            _context,                                       \
            _self.xyPlane,                                  \
            lineMesh,                                       \
            { 'marker'  : nasion,                           \
              'compare' : OBJECT_OT_SetNasion.fCompare }    \
        )

        # find the inion - too hard to distinguish
        '''
            # find the lowest vertex, as we did not remove vertices on
            # the back-side of the head the lowest vertex should be there
            # at the other line-beginning
        bpy.ops.object.mode_set(mode='EDIT')
        lineMesh = bmesh.from_edit_mesh( _self.xyPlane.data )

            # need to call this before direct indexing is possible
        if hasattr(lineMesh.verts, "ensure_lookup_table"):
            lineMesh.verts.ensure_lookup_table()

        lowestVertex = lineMesh.verts[0]
        for lineVertex in lineMesh.verts:
            lineVertex.select = False
            if lineVertex.co.z < lowestVertex.co.z:
                lowestVertex = lineVertex

            #climb upwards

        lowestVertex.select = True
        '''

        return {"FINISHED"}

class HelperGrid1020( object ):
    EEG1020_LUT = {
        "C3"    : None,
        "C4"    : None,
        "Cz"    : None,

        "F3"    : None,
        "F4"    : None,
        "F7"    : None,
        "F8"    : None,

        "Fpz"   : None,
        "Fp1"   : None,
        "Fp2"   : None,
        "Fz"    : None,

        "O1"    : None,
        "O2"    : None,
        "Oz"    : None,

        "P3"    : None,
        "P4"    : None,
        "Pz"    : None,

        "T3"    : None,
        "T4"    : None,
        "T5"    : None,
        "T6"    : None
    }

    def __init__( _self ):
        Utils.calcLineStats( bpy.context.scene.MyObjects[ OBJECT_OT_CreateNasionInionPlane.PLANE_NAME ] )
        Utils.calcLineStats( bpy.context.scene.MyObjects[ OBJECT_OT_CreateLRPlane.PLANE_NAME ] )

        print( "Initializing 10/20 helper grid..." )
        modelStats = bpy.context.scene.MyObjects[ "skin_stats" ]
        model      = bpy.context.scene.MyObjects[ "skin" ]


        # calculate alle the positions along the NI-line
        bpy.ops.object.mode_set(mode='OBJECT')
        niLine = bpy.context.scene.MyObjects[ OBJECT_OT_CreateNasionInionPlane.PLANE_NAME ]
        Utils.resetSelectedAndActiveObjects( )
        Utils.selectAndActivateObject( niLine ) 
        bpy.ops.object.mode_set(mode='EDIT')
        niLineMesh  = bmesh.from_edit_mesh( niLine.data )

        localCenter = niLine.matrix_world.inverted() * modelStats['center']                 # take axis inversion into consideration, see comment in OBJECT_OT_SetNasion for further details
        nasion = Utils.findLowestVertex( niLineMesh, lambda vertex: (localCenter.y > vertex.y) == bpy.context.scene.AxisInverted)

        _self.EEG1020_LUT[ "Fpz" ] =     \
            niLine.matrix_world * Utils.vertexAtLinestripFraction( niLineMesh, nasion, 0.1, bpy.types.Scene.MyObjects[ 'line_lengths' ][ niLine.name ] ).co

        _self.EEG1020_LUT[ "Fz" ] =     \
            niLine.matrix_world * Utils.vertexAtLinestripFraction( niLineMesh, nasion, 0.3, bpy.types.Scene.MyObjects[ 'line_lengths' ][ niLine.name ] ).co

        _self.EEG1020_LUT[ "Cz" ] =     \
            niLine.matrix_world * Utils.vertexAtLinestripFraction( niLineMesh, nasion, 0.5, bpy.types.Scene.MyObjects[ 'line_lengths' ][ niLine.name ] ).co

        _self.EEG1020_LUT[ "Pz" ] =     \
            niLine.matrix_world * Utils.vertexAtLinestripFraction( niLineMesh, nasion, 0.7, bpy.types.Scene.MyObjects[ 'line_lengths' ][ niLine.name ] ).co

        _self.EEG1020_LUT[ "Oz" ] =     \
            niLine.matrix_world * Utils.vertexAtLinestripFraction( niLineMesh, nasion, 0.9, bpy.types.Scene.MyObjects[ 'line_lengths' ][ niLine.name ] ).co

        bpy.ops.object.mode_set(mode='OBJECT')

        # calculate all the positions along the LR-line
        bpy.ops.object.mode_set(mode='OBJECT')
        lrLine = bpy.context.scene.MyObjects[ OBJECT_OT_CreateLRPlane.PLANE_NAME ]
        Utils.resetSelectedAndActiveObjects( )
        Utils.selectAndActivateObject( lrLine ) 
        bpy.ops.object.mode_set(mode='EDIT')
        lrLineMesh = bmesh.from_edit_mesh( lrLine.data )

        localCenter = lrLine.matrix_world.inverted() * modelStats['center']
        lPoint = Utils.findLowestVertex( lrLineMesh, lambda vertex: localCenter.x > vertex.x )

        _self.EEG1020_LUT[ "T3" ] =     \
            lrLine.matrix_world * Utils.vertexAtLinestripFraction( lrLineMesh, lPoint, 0.1, bpy.types.Scene.MyObjects[ 'line_lengths' ][ lrLine.name ] ).co

        _self.EEG1020_LUT[ "C3" ] =     \
            lrLine.matrix_world * Utils.vertexAtLinestripFraction( lrLineMesh, lPoint, 0.3, bpy.types.Scene.MyObjects[ 'line_lengths' ][ lrLine.name ] ).co

        _self.EEG1020_LUT[ "C4" ] =     \
            lrLine.matrix_world * Utils.vertexAtLinestripFraction( lrLineMesh, lPoint, 0.7, bpy.types.Scene.MyObjects[ 'line_lengths' ][ lrLine.name ] ).co

        _self.EEG1020_LUT[ "T4" ] =     \
            lrLine.matrix_world * Utils.vertexAtLinestripFraction( lrLineMesh, lPoint, 0.9, bpy.types.Scene.MyObjects[ 'line_lengths' ][ lrLine.name ] ).co


        # create the lower ring around the head
        horLineName = "horizontalHelperPlane"

        vec = _self.EEG1020_LUT[ "Fpz" ] - _self.EEG1020_LUT[ "Oz" ]
        ctr = _self.EEG1020_LUT[ "Oz" ] + 0.5 * vec
        ctr.x += modelStats[ 'dimensions' ].x * 0.5

        Utils.createPlaneBetweenPoints(
            [ _self.EEG1020_LUT[ "Fpz" ], _self.EEG1020_LUT[ "Oz" ] ],
            horLineName,
            ctr )

        horLine = bpy.context.scene.MyObjects[ horLineName ]
        Utils.resetSelectedAndActiveObjects( )
        Utils.selectAndActivateObject( horLine )
        bpy.ops.object.mode_set(mode='EDIT')
        horLineMesh = bmesh.from_edit_mesh( horLine.data )

        # another way to remove individual parts of the mesh
        # using the other way blender always crahsed, not using this way
        for face in horLineMesh.faces:
           face.select = True
        bpy.ops.mesh.delete(type='ONLY_FACE') #

        Utils.pojectFrameOntoSurface( model, horLine )

        # measure the positions along the just created ring
        Utils.calcLineStats( horLine )
        closest = Utils.findClosestMeshPoint( horLine, _self.EEG1020_LUT[ "Fpz" ] )
        horLineMesh = closest['mesh']

        _self.EEG1020_LUT[ "Fp2" ] =     \
            horLine.matrix_world * Utils.vertexAtLinestripFraction( horLineMesh, closest['point' ], 0.05, bpy.types.Scene.MyObjects[ 'line_lengths' ][ horLine.name ], 'x', True ).co

        _self.EEG1020_LUT[ "F8" ] =     \
            horLine.matrix_world * Utils.vertexAtLinestripFraction( horLineMesh, closest['point' ], 0.15, bpy.types.Scene.MyObjects[ 'line_lengths' ][ horLine.name ], 'x', True ).co

        _self.EEG1020_LUT[ "T6" ] =     \
            horLine.matrix_world * Utils.vertexAtLinestripFraction( horLineMesh, closest['point' ], 0.35, bpy.types.Scene.MyObjects[ 'line_lengths' ][ horLine.name ], 'x', True ).co

        _self.EEG1020_LUT[ "O2" ] =     \
            horLine.matrix_world * Utils.vertexAtLinestripFraction( horLineMesh, closest['point' ], 0.45, bpy.types.Scene.MyObjects[ 'line_lengths' ][ horLine.name ], 'x', True ).co

        _self.EEG1020_LUT[ "Fp1" ] =     \
            horLine.matrix_world * Utils.vertexAtLinestripFraction( horLineMesh, closest['point' ], 0.05, bpy.types.Scene.MyObjects[ 'line_lengths' ][ horLine.name ], 'x', False ).co

        _self.EEG1020_LUT[ "F7" ] =     \
            horLine.matrix_world * Utils.vertexAtLinestripFraction( horLineMesh, closest['point' ], 0.15, bpy.types.Scene.MyObjects[ 'line_lengths' ][ horLine.name ], 'x', False ).co

        _self.EEG1020_LUT[ "T5" ] =     \
            horLine.matrix_world * Utils.vertexAtLinestripFraction( horLineMesh, closest['point' ], 0.35, bpy.types.Scene.MyObjects[ 'line_lengths' ][ horLine.name ], 'x', False ).co

        _self.EEG1020_LUT[ "O1" ] =     \
            horLine.matrix_world * Utils.vertexAtLinestripFraction( horLineMesh, closest['point' ], 0.45, bpy.types.Scene.MyObjects[ 'line_lengths' ][ horLine.name ], 'x', False ).co

        # calculate the coords of the last four positions
        center = _self.EEG1020_LUT[ "F7" ] - _self.EEG1020_LUT[ "Fz" ]
        _self.EEG1020_LUT[ "F3" ] =     \
            _self.EEG1020_LUT[ "Fz" ] + 0.5 * center

        center = _self.EEG1020_LUT[ "F8" ] - _self.EEG1020_LUT[ "Fz" ]
        _self.EEG1020_LUT[ "F4" ] =     \
            _self.EEG1020_LUT[ "Fz" ] + 0.5 * center

        center = _self.EEG1020_LUT[ "T5" ] - _self.EEG1020_LUT[ "Pz" ]
        _self.EEG1020_LUT[ "P3" ] =     \
            _self.EEG1020_LUT[ "Pz" ] + 0.5 * center

        center = _self.EEG1020_LUT[ "T6" ] - _self.EEG1020_LUT[ "Pz" ]
        _self.EEG1020_LUT[ "P4" ] =     \
            _self.EEG1020_LUT[ "Pz" ] + 0.5 * center


    def getWorldCoord( _self, _1020coord ):
        try:
            return _self.EEG1020_LUT[ _1020coord ]
        except:
            print( "Invalid 10/20 coord provided!" )
            return None

class OBJECT_OT_ElectrodePrototype( bpy.types.Operator ):
    """Create the helper line to specify the nasion and the inion later"""
    mNumberRuns     = 0

    bl_idname       = "object.electrode_prototype"
    bl_label        = "Set Prototype"
    bl_options      = { "REGISTER", "UNDO" }
    bl_description  = "Erstellt und platziert eine prototypische Elektrode auf das Kopfmodell"

    def __init__( _self ):
        if not "1020grid" in bpy.context.scene.MyObjects:
            bpy.types.Scene.MyObjects["1020grid"] = HelperGrid1020()

        mNumberRuns = 0

    #
    # uses the helper grid to get the world coordinates of the 10/20 coordinate
    #
    #def obtainElectrodeCoord_from1020( _self, _1020coord ):
    #    coord = bpy.context.scene.MyObjects['1020grid'].getWorldCoord( _1020coord )
    #
    #    if coord:
    #        return coord
    #    else:
    #        return Vector( ( 23.0, 80.0, 220.0 ) )

    #
    # Assumptions on the state:
    # - helper grid exists
    #
    def execute( _self, _context ):
        skinSurface     = _context.scene.MyObjects["skin"]

        electrodePos    = Utils.obtainElectrodeCoord_from1020( _context.scene.ElectrodePos )
        electrodeX      = _context.scene.ElectrodeX
        electrodeY      = _context.scene.ElectrodeY
        ratio           = electrodeX / electrodeY
        electrodeName   = "electrode" + repr( _context.scene.IDCurrentElectrode )
        cuttingCubeName = electrodeName + "_cuttingcube"

        meshStat  = _context.scene.MyObjects["skin_stats"]
        origin    = meshStat['origin']
        dim       = meshStat['dimensions']
        center    = meshStat['center']
        normalVec = electrodePos - center
        direction = Vector( -1 * normalVec )
        direction.normalize()
        normalVec.normalize()

        if electrodePos == None:
            bpy.ops.error.message('INVOKE_DEFAULT',
                                  mType = "Error",
                                  mMessage = "Invalid electrode position provided!")
            return { "FINISHED" }

        # set skin surface to object mode
        bpy.ops.object.mode_set(mode='OBJECT')
        
        smoothSkinSurface  = _context.scene.MyObjects["smoothskin"]
        
        # now obtain the closest normal to the electrode pos
        rayStart = electrodePos + normalVec*100
        (did_hit,hit_loc,hit_face_normal,face_index) \
            = smoothSkinSurface.ray_cast(rayStart, direction)
        
        # check if the provided electrode-ID already exists
        # yes : delete the existing electrode with this ID, create a new one
        # no  : create a new electrode and increment the IDs accordingly
        if( bpy.data.objects.get( cuttingCubeName ) is not None ):
            Utils.resetSelectedAndActiveObjects()
            Utils.selectAndActivateObject( _context.scene.MyObjects[ cuttingCubeName ] )
            bpy.ops.object.delete()

        Utils.resetSelectedAndActiveObjects()
        ############### State changed #################
        # ACTIVE: None, SELECTED : None, MODE: Object #
        ###############################################

        # create a helper vertex, at the position of the electrode
        # we'll use this point to correctly orient the cutting-cube before
        # projecting it onto the skin surface
        helper = _context.scene.MyObjects["helper"]
        Utils.resetSelectedAndActiveObjects()
        Utils.selectAndActivateObject( helper )
        ############### State changed ######################
        # ACTIVE: helper, SELECTED : helper, MODE: Object  #
        ####################################################

        # create the cube use to cut the electrode contact-surface from the surface
        bpy.ops.mesh.primitive_cube_add( radius=electrodeX/2, location=electrodePos )
        cuttingCube       = bpy.context.active_object
        cuttingCube.name  = cuttingCubeName

        ###################### State changed ########################
        # ACTIVE: cuttingcube, SELECTED : cuttingcube, MODE: Object #
        #############################################################

        # orient cube to helper object (according to the normal of the closest face)
        #   -> in the standard-case electrode positions
        #   -> positions that are not along the center-z-axis
        quat_angle = hit_face_normal.to_track_quat('Y','Z')
            # due to a restriction in blender concerning mesh updates
            # (see: http://blender.stackexchange.com/a/44324)
            # we cannot do this:
            #   cuttingCube.dimensions[0] = electrodeX
            #   cuttingCube.dimensions[electrodeYDim] = electrodeY
            #   cuttingCube.dimensions[electrodeZDim] = min( min( electrodeX * 1.5, electrodeY * 1.5), 100 )
            # but we have to assign the values at once
        cuttingCube.dimensions = electrodeX, min( min( electrodeX * 1.5, electrodeY * 1.5), 100 ), electrodeY

            # special case for positions along the center-z-axis
        center_coordinate = re.compile("^.*z$")
        if center_coordinate.match(_context.scene.ElectrodePos):    
            # the track axis changes for the very center position and therefor also the 
            # 'depth-axis' of the cube (the axis directed towards the head)
            if _context.scene.ElectrodePos == "Cz":
                quat_angle = hit_face_normal.to_track_quat('Z','Y')
                cuttingCube.dimensions = electrodeX, electrodeY, min( min( electrodeX * 1.5, electrodeY * 1.5), 100 ) 
            # for all the positions along the center-z-axis we have to reset the z-rotation, because
            # strangley when tracking the surface-normal of the head, the cube always gets slightly
            # rotated around the z-axis as well (I suspect a kind of 'gimbal lock' here...)
            euler_angle = quat_angle.to_euler()
            euler_angle.z = 0.0
            quat_angle = euler_angle.to_quaternion()
            
        bpy.context.object.rotation_mode = 'QUATERNION'
        bpy.context.object.rotation_quaternion = quat_angle

            # project cube onto surface (without warping it, note that this
        # is a 'constraint' not a modifier)
        bpy.ops.object.constraint_add( type='SHRINKWRAP' )
        swc        = cuttingCube.constraints[0]
        swc.target = skinSurface
        swc.shrinkwrap_type = "NEAREST_SURFACE"    #to make the cube stick to the surface
        swc.project_axis    = 'POS_Y'


        _context.scene.MyObjects[ cuttingCubeName ] = cuttingCube
        
        localCenter = cuttingCube.matrix_world.inverted() * cuttingCube.location
                
        # keep this for later reference:
        #   here we can assign an individual color to each of the vertices of the mesh
        '''
        cuttingCubeMesh = cuttingCube.data 
        
        if cuttingCubeMesh.vertex_colors:
            vcolLayer = cuttingCubeMesh.vertex_colors.active
        else:
            vcolLayer = cuttingCubeMesh.vertex_colors.new()
            
        for poly in cuttingCubeMesh.polygons:       # all polygons associated with this mesh-object --> faces (?)
            for loopIndex in poly.loop_indices:     # all indices (of vertices) of the current polygon-corner (=closed loop of edges); note that we have in total 24 vertizes for a cube, because each 
                loopVertexIndex = cuttingCubeMesh.loops[loopIndex].vertex_index
                vertex = cuttingCubeMesh.vertices[loopVertexIndex]
                localCenter = cuttingCube.matrix_world.inverted() * cuttingCube.location
                if vertex.co.x < localCenter.x:                         # @TODO adapt for Cz, in this case it mus be the y-coordinate              
                    vcolLayer.data[loopIndex].color = (1.0, 0.0, 1.0)
                    print("painting vertex ",loopIndex)
                else:
                    print( vcolLayer.data[loopIndex].color )
        
        
            # this is NOT enough: we have to address each vertex for each face individually
            #for vertex in cuttingCubeMesh.vertices:
            #    color = vertex.normal
            #    vcolLayer.data[vertex.index].color = color
            #    print("painting vertex ",vertex.index," to color ", color[0], color[1], color[2])
            #
        
        cuttingCubeMesh.vertex_colors.active = vcolLayer
        cuttingCubeMesh.update()
        
        # enable rendering the vertex paint in object-mode (not only in mode: 'VERTEX_PAINT', bpy.ops.object.mode_set(mode='VERTEX_PAINT');)
        bpy.context.space_data.show_textured_solid = True
        '''
        
        #~ if (_context.scene.AutoRun ):
            #~ if( len( _context.scene.PredefinedElectrodeProperties ) > _context.scene.IDCurrentElectrode ):
                #~ _context.scene.ElectrodePos = _context.scene.PredefinedElectrodeProperties[ _context.scene.IDCurrentElectrode ]['coord']
                #~ _context.scene.ElectrodePos = _context.scene.PredefinedElectrodeProperties[ _context.scene.IDCurrentElectrode ]['coord']
                #~ _context.scene.ElectrodeX   = _context.scene.PredefinedElectrodeProperties[ _context.scene.IDCurrentElectrode ]['size'][0]
                #~ _context.scene.ElectrodeY   = _context.scene.PredefinedElectrodeProperties[ _context.scene.IDCurrentElectrode ]['size'][1]

        return {"FINISHED"}

#
# main operator that creates the electrode-surfaces
#
class OBJECT_OT_MountElectrodes( bpy.types.Operator ):
    """Mount electrodes onto skin surface"""

    bl_idname       = "object.mount_electrodes"
    bl_label        = "Mount Electrodes"
    bl_options      = { "REGISTER", "UNDO" }
    bl_description  = "Setzt die Elektroden auf den ausgewählten Punkt"

    mElectrodeParts = [ "contact", "walls" ]
    mNumberRuns     = 0

    DEFAULT_ELECTRODE_SIZE_X   = 50
    DEFAULT_ELECTRODE_SIZE_Y   = 50
    DEFAULT_ELECTRODE_HEIGHT   = 5
    DEFAULT_ELECTRODE_POSITION = "Cz"

    def __init__ ( _self ):
        if not "1020grid" in bpy.context.scene.MyObjects:
            bpy.types.Scene.MyObjects["1020grid"] = HelperGrid1020()

    def cleanExistingElectrode( _self, _context, _electrodeName, _wallsSurfName, _contactSurfName, _gelLayerName ):
        Utils.resetSelectedAndActiveObjects()
        Utils.selectAndActivateObject( _context.scene.MyObjects[ _contactSurfName  ] )
        bpy.ops.object.delete()
        Utils.resetSelectedAndActiveObjects()
        Utils.selectAndActivateObject( _context.scene.MyObjects[ _electrodeName  ] )
        bpy.ops.object.delete()
        Utils.resetSelectedAndActiveObjects()
        Utils.selectAndActivateObject( _context.scene.MyObjects[ _wallsSurfName ] )
        bpy.ops.object.delete()
        # we may or may not have a gel layer created previously - so check first!
        if( bpy.data.objects.get(_gelLayerName) is not None ):
            Utils.resetSelectedAndActiveObjects()
            Utils.selectAndActivateObject( _context.scene.MyObjects[ _gelLayerName ] )
            bpy.ops.object.delete()



    def execute( _self, _context ):
        if bpy.context.scene.UseSmootSkin:
            skinSurface     = _context.scene.MyObjects["smoothskin"]
        else:
            skinSurface     = _context.scene.MyObjects["skin"]

        electrodeName   = "electrode" + repr( _context.scene.IDCurrentElectrode )
        cuttingCubeName = electrodeName + "_cuttingcube"
        contactSurfName = electrodeName + "_contact"
        wallsSurfName   = electrodeName + "_walls"
        gelLayerName    = electrodeName + '_gel'
        
        cuttingCube = _context.scene.MyObjects[ cuttingCubeName ]
        
        electrodePos = Utils.obtainElectrodeCoord_center( cuttingCube )
        meshStat     = _context.scene.MyObjects["skin_stats"]
        dim          = meshStat['dimensions']
        center       = meshStat['center']
        normalVec    = electrodePos - center
        direction    = Vector( -1 * normalVec )
        direction.normalize()
        normalVec.normalize()

        if( bpy.data.objects.get(electrodeName) is not None ):
            _self.cleanExistingElectrode( _context, electrodeName, wallsSurfName, contactSurfName, gelLayerName )

        ### 1: prepare skin surface for CSG by duplicating it.
        # This way we retain the original surface and cut a copy
        Utils.resetSelectedAndActiveObjects()
        Utils.selectAndActivateObject( skinSurface )

        # ! STATE CHANGED: ACTIVE: skin, SELECTED : skin, MODE: Object

        bpy.ops.object.duplicate()
        electrodeContactSurface = bpy.context.active_object
        electrodeContactSurface.name = contactSurfName

        # ! STATE CHANGED:  ACTIVE: contact, SELECTED : contact, MODE: Object

        ### 2: cut the cuttingcube out of the just duplicated skin-surface
        # The resulting intersection will be the contact surface of the electrode
        bpy.ops.object.modifier_add( type='BOOLEAN' )
        csg             = electrodeContactSurface.modifiers[0]
        csg.object      = cuttingCube
        csg.operation   = "INTERSECT"

            # in order to apply a modifier two preconditions have to satisfied:
            # a) the object whose modifiers should be applied needs to be the active object
            # b) the name of the modifier to be applied has to be provided (e.g. csg.name)
        bpy.ops.object.modifier_apply( apply_as='DATA', modifier=csg.name )


            # remove all faces except the outer part of the cut,
            # as this is the real contact surface
        bpy.ops.object.mode_set(mode='EDIT')

        # ! STATE CHANGED: ACTIVE: contact, SELECTED : contact, MODE: Edit

            # remove any faces whose normals point inside the skin
        Utils.removeVerticesAccordingToDirection( electrodeContactSurface,  0.1, direction )

            # remove any faces that are not connected to the frontmost vertex
            # (removes 'noise' due to internal geometry intersected by the cube)
        electrodeContactSurfaceMesh = bmesh.from_edit_mesh( electrodeContactSurface.data )
        outmostVert = None
        largestDist = -1 * (sys.maxsize) - 1
        for vertex in electrodeContactSurfaceMesh.verts:
            dist = Utils.euklidianDist( center, vertex.co )
            if( dist > largestDist ):
                outmostVert = vertex
                largestDist = dist

        outmostVert.select = True

        bpy.ops.mesh.select_linked()
        bpy.ops.mesh.select_all(action='INVERT')

        bpy.ops.mesh.delete(type='VERT')                  # True = make changes to the mesh persitent, recalculate normals
        bmesh.update_edit_mesh(electrodeContactSurface.data, True)

            # select the remaining vertices again, fill any holes that might
            # be created due to skew faces
        bpy.ops.mesh.select_all(action='TOGGLE')
        bpy.ops.mesh.fill_holes(sides=16)

        bpy.ops.object.mode_set(mode='OBJECT')
        bpy.ops.object.origin_set( type='ORIGIN_GEOMETRY' )

        # ! STATE CHANGED: ACTIVE: contact, SELECTED : contact, MODE: Object

        _context.scene.MyObjects[ contactSurfName ] = electrodeContactSurface

        ### 3. if desired create the gel layer
        # this is done the same way, the electrode will be created
        # as an extrusion of the contact-surface
        if bpy.context.scene.UseGelLayer:
            Utils.resetSelectedAndActiveObjects()
            Utils.selectAndActivateObject( electrodeContactSurface )
            bpy.ops.object.duplicate()

            gelLayer = bpy.context.active_object
            gelLayer.name = gelLayerName

            Utils.resetSelectedAndActiveObjects()
            Utils.selectAndActivateObject( gelLayer )
            bpy.ops.object.mode_set(mode='EDIT')

            # ! STATE CHANGED: ACTIVE: gel, SELECTED : gel, MODE: Edit
            extrusionVector = (0.5*normalVec.x, 0.5*normalVec.y, 0.5*normalVec.z)

            bpy.ops.mesh.extrude_region_move(
                TRANSFORM_OT_translate={"value":extrusionVector}
            )
            bpy.ops.object.mode_set(mode='OBJECT')
            bpy.ops.object.origin_set( type='ORIGIN_GEOMETRY' )


            # ! STATE CHANGED: ACTIVE: gel, SELECTED : gel, MODE: Object
            _context.scene.MyObjects[ gelLayerName ] = gelLayer

            Utils.resetSelectedAndActiveObjects()
            Utils.selectAndActivateObject( electrodeContactSurface )
            electrodeContactSurface.location[0] += extrusionVector[0]
            electrodeContactSurface.location[1] += extrusionVector[1]
            electrodeContactSurface.location[2] += extrusionVector[2]



        ### 4. now create the full electrode from the contact surface
        # This is done by extruding the contact surface
            # extrude the mesh to get the full electrode,
            # again we work on a copy here,
            # as we also want to keep the 'contact_surface'
        Utils.resetSelectedAndActiveObjects()
        Utils.selectAndActivateObject( electrodeContactSurface )
        bpy.ops.object.duplicate()
        electrodeFull = bpy.context.active_object
        electrodeFull.name = electrodeName

        Utils.resetSelectedAndActiveObjects()
        Utils.selectAndActivateObject( electrodeFull )
        bpy.ops.object.mode_set(mode='EDIT')

        # ! STATE CHANGED: ACTIVE: full, SELECTED : full, MODE: Edit

        for i in range( 0, int( bpy.context.scene.ElectrodeZ) ):
            bpy.ops.mesh.extrude_region_move(
                TRANSFORM_OT_translate={"value":(normalVec.x, normalVec.y, normalVec.z)}
            )

        bpy.ops.object.mode_set(mode='OBJECT')
        bpy.ops.object.origin_set( type='ORIGIN_GEOMETRY' )


        # ! STATE CHANGED: ACTIVE: full, SELECTED : full, MODE: Object

        _context.scene.MyObjects[ electrodeName ] = electrodeFull

        #### 5. now create the electrode 'walls' as a separate mesh
            # remove contact_surface faces from 'full_electrode' to get the
            # walls faces -> angain copy the full electrode first
        bpy.ops.object.duplicate()
        electrodeWalls = bpy.context.active_object
        electrodeWalls.name = wallsSurfName

        Utils.resetSelectedAndActiveObjects()
        Utils.selectAndActivateObject( electrodeWalls )
        bpy.ops.object.mode_set(mode='EDIT')

        # ! STATE CHANGED: ACTIVE: walls, SELECTED : walls, MODE: Edit

        Utils.removeFacesAccordingToDirection( electrodeWalls, 0.25, direction )

        bpy.ops.object.mode_set(mode='OBJECT')

        # ! STATE CHANGED: ACTIVE: walls, SELECTED : walls, MODE: Object

        _context.scene.MyObjects[ wallsSurfName ] = electrodeWalls
        bpy.ops.object.origin_set( type='ORIGIN_GEOMETRY' )

        _context.scene.IDCurrentElectrode += 1
 
        # populate input fields with default values OR values from input string for the next defined electrode (if present)       
        Utils.populateInputs( _context, _context.scene.IDCurrentElectrode )

        return {"FINISHED"}

if __name__ == "__main__":
    unregister()
    register()

