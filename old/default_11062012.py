import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import savefig
#from matplotlib.figure import Figure
#from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import networkx as nx
import numpy as np
import scipy
#import scipy.stats
exec('from applications.%s.modules import core as core'%request.application)
exec('from applications.%s.modules import community as community'%request.application)
#exec('from applications.%s.modules import plot_network as plot_network'%request.application)
#exec('from applications.%s.controllers.appadmin import *'%request.application)
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch, cm
from reportlab.lib import colors, utils
from reportlab.lib.pagesizes import letter
import re
import time
import datetime
from gluon.storage import Storage

def user():
    return dict(form=auth())

def study_networks():
    if auth.user:
        rows_allowed = db((db.upload_data.share=='Public') | (db.upload_data.email==auth.user.email)).select()
    else:
        rows_allowed = db((db.upload_data.share=='Public')).select()
    study_names = list(set([row.study_name for row in rows_allowed]))
    study_networks = ''
    if request.vars.study_name or request.vars.study_name_1 or request.vars.study_name_2:
        for row in rows_allowed:
            if (row.study_name==request.vars.study_name) or (row.study_name==request.vars.study_name_1) or (row.study_name==request.vars.study_name_2):
                study_networks = study_networks + '<option value="%s">%s</option>' %(row.network_name, row.network_name)
    return study_networks

def webgl_html():
    network = db(db.upload_data.network_name==request.vars.network_name).select().first()
    matfilename = network.connectivity_matrix_file
    matfilepath = os.path.join(request.folder,'uploads',matfilename)
    centersfilepath = os.path.join(request.folder,'uploads',network.region_xyz_centers_file)
    #molfile_str = plot_matrix_sim_molfile(matfilepath, centersfilepath, int(request.vars.chosen_density))
    weight = request.vars.weight
    weight_edges = (weight!='Binary')
    if weight == 'Binary':
        binarize=True
    else:
        binarize=False
    chosen_density = int(request.vars.chosen_density)
    node_metric = request.vars.node_metric
    webgl_url, node_text, edge_text = plot_matrix_webgl(matfilepath,centersfilepath, network, weight=weight, chosen_density=chosen_density, node_metric=node_metric)
    return dict(webgl_url=webgl_url, node_text=node_text, edge_text=edge_text)

def plot_matrix_webgl(connectmat_file,centers_file,network, chosen_density, node_metric='deg',
                      weight='Binary',node_scale_factor=.25, edge_radius=.05, names_file=None):
#def plot_matrix_webgl(weight='Binary',
#                      weight_edges=0, node_scale_factor=.25, edge_radius=.05, names_file=None):
    """Given a connectivity matrix and a (x,y,z) centers file for each region,
    generate a chemdoodle webgl 3d javascript"""
    
    #######
    #network = db(db.upload_data.network_name==request.vars.network_name).select().first()
    #matfilename = network.connectivity_matrix_file
    #matfilepath = os.path.join(request.folder,'uploads',matfilename)
    #centersfilepath = os.path.join(request.folder,'uploads',network.region_xyz_centers_file)
    ##molfile_str = plot_matrix_sim_molfile(matfilepath, centersfilepath, int(request.vars.chosen_density))
    #weight = request.vars.weight
    weight_edges = (weight!='Binary')
    if weight == 'Binary':
        binarize=True
    else:
        binarize=False
    #connectmat_file = matfilepath
    #centers_file = centersfilepath
    #if int(request.vars.chosen_density) > 20:
    if chosen_density > 20:
        threshold_pct = 20
    else:
        #threshold_pct = int(request.vars.chosen_density)
        threshold_pct = int(chosen_density)
    ########
    
    matrix = core.file_reader(connectmat_file)
    matrix_array = np.array(matrix)
    nodes = core.file_reader(centers_file)
    #metric = request.vars.node_metric
    
    webgl_filename = os.path.join(network.network_name \
                     + '_' + weight + '_' + str(threshold_pct) \
                     + '_' + node_metric + '_' + 'webgl' + '.js')
    webgl_filepath = os.path.join(request.folder, 'static', webgl_filename)
    
    #if os.path.exists(webgl_filepath):
    #    print '1'
    #    webgl_url = URL(request.application, 'static', webgl_filename)
    #    return webgl_url
    
    node_text = ''
    if node_metric:
        if node_metric == 'deg':
            metric_num = 1
            node_color = "orange"
            node_text = 'Node radius is currently based on the node\'s <font color="orange">degree</b></font>.'
        if node_metric == 'bc':
            metric_num = 2
            node_color = "chartreuse"
            node_text = 'Node radius is currently based on the node\'s <font color="#7fff00">betweenness centrality</b></font>.'
        if node_metric == 'cc':
            metric_num = 3
            node_color = "aqua"
            node_text = 'Node radius is currently based on the node\'s <font color="#00ffff"><b>clustering coefficient</b></font>.'
        if node_metric == 'mod':
            metric_num = 4
            node_text = 'Node color is currently based on the node\'s <b>module membership</b>.'
    else:
        metric_num = 0
    
    f = open(webgl_filepath, 'w')
    f.write("var transform3d = new ChemDoodle.TransformCanvas3D('transformBallAndStick', 500, 500);\n")
    f.write("transform3d.specs.set3DRepresentation('Ball and Stick');\n")
    f.write("transform3d.specs.backgroundColor = 'white';\n")
    f.write("ChemDoodle.ELEMENT['C'].jmolColor = '#99FF00';\n")
    f.write("var molecule = new ChemDoodle.structures.Molecule()\n")
    
    if names_file:
        names = core.file_reader(names_file,1)
    num_nodes = len(nodes)
    
    edge_thresh_pct = threshold_pct / 100.0
    edge_text = 'Displaying the top <b>%s%%</b> strongest connections' %threshold_pct
    matrix_flat=np.array(matrix).flatten()
    matrix_flat[np.isinf(matrix_flat)] = 0
    matrix_flat[np.isnan(matrix_flat)] = 0
    edge_thresh=np.sort(matrix_flat)[len(matrix_flat)-int(len(matrix_flat)*edge_thresh_pct)]
    num_edges = len(matrix_array[np.nonzero(matrix_array > edge_thresh)])/2
    cutoff = 800
    if num_edges > cutoff:
        edge_thresh = np.sort(matrix_flat)[len(matrix_flat) - cutoff]
        num_edges = cutoff
        edge_thresh_pct_cap = int(((num_edges*2)/float(len(matrix_flat)))*100)
        edge_text = 'Displaying the top <b>%s%%</b> strongest connections (capped to enable rendering).' %edge_thresh_pct_cap
    
    # maximum spanning tree
    inv_mat = np.matrix(matrix)
    inv_mat[np.nonzero(inv_mat)] = np.max(inv_mat)+.01-inv_mat[np.nonzero(inv_mat)]
    Ginv = nx.from_numpy_matrix(inv_mat)
    mst = nx.minimum_spanning_edges(Ginv)
    edgelist = list(mst)
    edgelist_ids = [(x[0],x[1]) for x in edgelist]
    show_mst = False
    #####
    
    edge_min = min(matrix_flat)
    edge_max = max(matrix_flat)
    edge_range = edge_max - edge_min
    
    nodes_min = min(np.array(nodes).flatten())
    nodes = nodes + (-nodes_min)
    nodes_max = max(nodes.flatten())
    nm10 = nodes_max / 100
    nodes = nodes / (nm10 + .01)
    
    edge_scale_factor = .09 # .19
    node_scale_factor = .17 # defaults are .35
    
    metrics = session.metrics
    if metric_num:
        #if metrics[0] == ['Region Name', 'Degree', 'Clustering Coefficient', 'Betweenness Centrality', 'Module']:
        #    metrics = metrics[1:]
        metric_all = [x[metric_num] for x in metrics]
        metric_max = max(metric_all)
        # scale node radius to between .1 and .25
        metric_all_scaled = [(((x/float(metric_max))*node_scale_factor)+.05) for x in metric_all]

    for count,(x,y,z) in enumerate(nodes):
        # nodes are named c0, c1, c2, c3...
        if metric_num:
            if node_metric == 'mod':
                f.write("var c%s = new ChemDoodle.structures.Atom('C',%s,%s,%s,%s,'%s');\n" %(str(count), x/10, y/10, z/10, node_scale_factor, session.module_colors[count]))
            else:
                f.write("var c%s = new ChemDoodle.structures.Atom('C',%s,%s,%s,%s,'%s');\n" %(str(count), x/10, y/10, z/10, metric_all_scaled[count], node_color))
        else:
            f.write("var c%s = new ChemDoodle.structures.Atom('C',%s,%s,%s,%s);\n" %(str(count), x/10, y/10, z/10, node_scale_factor))
        f.write("molecule.atoms[%s]=c%s;\n" %(count, str(count)))
        if names_file:
            pass
    count = 0
    for i in range(num_nodes-1):
        x0,y0,z0=nodes[i]
        for j in range(i+1, num_nodes):
            if show_mst:
                if (i,j) in edgelist_ids:
                    x1,y1,z1=nodes[j]
                    if weight_edges:
                        node_weight = ((matrix[i][j] / edge_max) * edge_scale_factor) + .01
                        f.write("var b%s = new ChemDoodle.structures.Bond(c%s, c%s,1, %s);\n" %(str(count), str(i), str(j), node_weight))
                    else:
                        f.write("var b%s = new ChemDoodle.structures.Bond(c%s, c%s,1, %s);\n" %(str(count), str(i), str(j), edge_radius))
                    f.write("molecule.bonds[%s]=b%s;\n" %(count, str(count)))
                    count += 1
            else:
                if matrix[i][j] > edge_thresh:
                    x1,y1,z1=nodes[j]
                    if weight_edges:
                        node_weight = ((matrix[i][j] / edge_max) * edge_scale_factor) + .01
                        f.write("var b%s = new ChemDoodle.structures.Bond(c%s, c%s,1, %s);\n" %(str(count), str(i), str(j), node_weight))
                    else:
                        f.write("var b%s = new ChemDoodle.structures.Bond(c%s, c%s,1, %s);\n" %(str(count), str(i), str(j), edge_radius))
                    f.write("molecule.bonds[%s]=b%s;\n" %(count, str(count)))
                    count += 1
    f.write("transform3d.loadMolecule(molecule);\n")
    f.close()
    webgl_url = URL(request.application, 'static', webgl_filename)
    #webgl_html = "<script type='text/javascript' src='%s'></script>" %(webgl_url)
    return webgl_url, node_text, edge_text

def trunc(f, n):
    '''Truncates/pads a float f to n decimal places without rounding'''
    return ('%.*f' % (n + 1, f))[:-1]

def index():
    if auth.user:
        rows_allowed = db((db.upload_data.share=='Public') | (db.upload_data.email==auth.user.email)).select()
        study_names = list(set([row.study_name for row in rows_allowed]))
        if request.vars.study_name:
            study_networks = [row.network_name for row in rows_allowed if row.study_name==request.vars.study_name]
        elif request.vars.study_name_cur:
            study_networks = [row.network_name for row in rows_allowed if row.study_name==request.vars.study_name_cur]
        else:
            study_networks = ''
        form = SQLFORM.factory(
            Field('study_name', requires=IS_IN_SET(study_names),
                  default = (lambda x: x if x else '')(request.vars.study_name_cur)),
            Field('network_name',
                  requires=IS_IN_SET(study_networks),
                  #db.upload_data,
                  #requires=IS_IN_DB(db((db.upload_data.share=='Public') | (db.upload_data.email==auth.user.email)),
                  #                  'upload_data.id',
                  #                  '%(network_name)s'),
                  default = (lambda x: x if x else '')(request.vars.network_name_cur)),
            Field('weight',requires=IS_IN_SET([('Binary','Binary'),('Weighted','Weighted')]), default='Binary'),
            Field('edge_density',default=20,requires=IS_INT_IN_RANGE(0, 101, error_message=T('Must be integer in range 0-100'))),
            Field('orientation',requires=IS_IN_SET([('Axial','Axial'),('Sagittal','Sagittal'),('Coronal','Coronal')]), default='Axial'),
            labels = {'weight':'Weighting scheme:',
                      'edge_density':'% of edges to include:'},
            submit_button='Analyze',
            keepopts=['network_name', 'weight', 'edge_density'])
    else:
        rows_allowed = db((db.upload_data.share=='Public')).select()
        study_names = list(set([row.study_name for row in rows_allowed]))
        if request.vars.study_name:
            study_networks = [row.network_name for row in rows_allowed if row.study_name==request.vars.study_name]
        elif request.vars.study_name_cur:
            study_networks = [row.network_name for row in rows_allowed if row.study_name==request.vars.study_name_cur]
        else:
            study_networks = ''
        form = SQLFORM.factory(
            Field('study_name', requires=IS_IN_SET(study_names),
                  default = (lambda x: x if x else '')(request.vars.study_name_cur)),
            Field('network_name',
                  db.upload_data,
                  requires=IS_IN_SET(study_networks),
                  #requires=IS_IN_DB(db(db.upload_data.share=='Public'),
                  #                  'upload_data.id',
                  #                  '%(network_name)s'),
                  default = (lambda x: x if x else '')(request.vars.network_name_cur)),
            Field('weight',requires=IS_IN_SET([('Binary','Binary'),('Weighted','Weighted')]), default='Binary'),
            Field('edge_density',default=20,requires=IS_INT_IN_RANGE(0, 101, error_message=T('Must be integer in range 0-100'))),
            Field('orientation',requires=IS_IN_SET([('Axial','Axial'),('Sagittal','Sagittal'),('Coronal','Coronal')]), default='Axial'),
            labels = {'weight':'Weighting scheme:',
                      'edge_density':'% of edges to include:'},
            submit_button='Analyze',
            keepopts=['network_name', 'weight', 'edge_density'])
    form[0][0].insert(-1,"Choose from the list of available studies")
    form[0][1].insert(-1,"Choose from the list of shared networks for the chosen study")
    form[0][2].insert(-1,"Choose a weighting scheme for the network, either weighted or binarized edges")
    form[0][3].insert(-1,"Choose the percent of connections to keep (0-100%)")
    form[0][4].insert(-1,"Choose the orientation from which to view the network")
    if form.accepts(request.vars, session):
        cache.ram.clear()
        redirect(URL(r=request,f='single_view',vars=form.vars))
    return dict(form=form)

def lesion():
    if auth.user:
        rows_allowed = db((db.upload_data.share=='Public') | (db.upload_data.email==auth.user.email)).select()
        study_names = list(set([row.study_name for row in rows_allowed]))
        if request.vars.study_name:
            study_networks = [row.network_name for row in rows_allowed if row.study_name==request.vars.study_name]
        elif request.vars.study_name_cur:
            study_networks = [row.network_name for row in rows_allowed if row.study_name==request.vars.study_name_cur]
        else:
            study_networks = ''
        form = SQLFORM.factory(
            Field('study_name', requires=IS_IN_SET(study_names),
                  default = (lambda x: x if x else '')(request.vars.study_name_cur)),
            Field('network_name',
                  db.upload_data,
                  requires=IS_IN_SET(study_networks),
                  default = (lambda x: x if x else '')(request.vars.network_name_cur)),
            Field('weight',requires=IS_IN_SET([('Binary','Binary'),('Weighted','Weighted')]), default='Binary'),
            Field('edge_density',default=20,requires=IS_INT_IN_RANGE(0, 101, error_message=T('Must be integer in range 0-100'))),
            Field('orientation',requires=IS_IN_SET([('Axial','Axial'),('Sagittal','Sagittal'),('Coronal','Coronal')]), default='Axial'),
            labels = {'weight':'Weighting scheme:',
                      'edge_density':'% of edges to include:'},
            submit_button='Get Regions',
            keepopts=['study_name', 'network_name', 'weight', 'edge_density', 'orientation'])
    else:
        rows_allowed = db((db.upload_data.share=='Public')).select()
        study_names = list(set([row.study_name for row in rows_allowed]))
        if request.vars.study_name:
            study_networks = [row.network_name for row in rows_allowed if row.study_name==request.vars.study_name]
        elif request.vars.study_name_cur:
            study_networks = [row.network_name for row in rows_allowed if row.study_name==request.vars.study_name_cur]
        else:
            study_networks = ''
        form = SQLFORM.factory(
            Field('study_name', requires=IS_IN_SET(study_names),
                  default = (lambda x: x if x else '')(request.vars.study_name_cur)),
            Field('network_name',
                  db.upload_data,
                  requires=IS_IN_SET(study_networks),
                  default = (lambda x: x if x else '')(request.vars.network_name_cur)),
            Field('weight',requires=IS_IN_SET([('Binary','Binary'),('Weighted','Weighted')]), default='Binary'),
            Field('edge_density',default=20,requires=IS_INT_IN_RANGE(0, 101, error_message=T('Must be integer in range 0-100'))),
            Field('orientation',requires=IS_IN_SET([('Axial','Axial'),('Sagittal','Sagittal'),('Coronal','Coronal')]), default='Axial'),
            labels = {'weight':'Weighting scheme:',
                      'edge_density':'% of edges to include:'},
            submit_button='Get Regions',
            keepopts=['study_name', 'network_name', 'weight', 'edge_density', 'orientation'])
    form[0][0].insert(-1,"Choose from the list of available studies")
    form[0][1].insert(-1,"Choose from the list of shared networks")
    form[0][2].insert(-1,"Choose a weighting scheme for the network, either weighted or binarized edges")
    form[0][3].insert(-1,"Choose the percent of connections to keep (0-100%)")
    form[0][4].insert(-1,"Choose the orientation from which to view the network")
    #form = SQLFORM.factory(
    #    Field('network_name',db.upload_data,requires=IS_IN_DB(db,'upload_data.id','%(network_name)s')),
    #    Field('weight',requires=IS_IN_SET([('Binary','Binary'),('Weighted','Weighted')]), default='Weighted'),
    #    Field('edge_density',default=100,requires=IS_INT_IN_RANGE(0, 101, error_message=T('Must be integer in range 0-100'))),
    #    Field('orientation',requires=IS_IN_SET([('Axial','Axial'),('Sagittal','Sagittal'),('Coronal','Coronal')]), default='Axial'),
    #    submit_button='Get Regions')
    #form[0][0].insert(-1,"Choose from the list of shared networks")
    #form[0][1].insert(-1,"Choose a weighting scheme for the network, either weighted or binarized edges")
    #form[0][2].insert(-1,"Choose the percent of connections to keep (0-100%)")
    fullnames_filepath = None
    form2 = FORM()
    if form.accepts(request.vars, session, keepvalues=True):
        session.lesion = {}
        session.lesion['network_1_name'] = request.vars.network_name
        session.lesion['network_2_name'] = request.vars.network_name
        session.lesion['weight_1'] = request.vars.weight
        session.lesion['weight_2'] = request.vars.weight
        session.lesion['edge_density_1'] = request.vars.edge_density
        session.lesion['edge_density_2'] = request.vars.edge_density
        session.lesion['orientation'] = request.vars.orientation
        
        network = db(db.upload_data.network_name==request.vars.network_name).select().first()
        fullnames_filepath = os.path.join(request.folder,'uploads',network.region_names_full_file)
        session.fullnames = file_reader(fullnames_filepath,1)
        
        components = [LI(INPUT(_type="checkbox", _name=region), region) for region in session.fullnames]
        components.append(INPUT(_type='submit', _value='Go')) #, _action=URL(r=request,f='compare_view',vars=vars)))
        form2 = FORM(*components)
    if form2.accepts(request.vars, session):
        lesion_regions = [y for y in request.vars.keys() if request.vars[y]=='on']
        session.lesion['region_nums'] = [session.fullnames.index(x) for x in lesion_regions] # region numbers
        cache.ram.clear()
        redirect(URL(f='compare_view', args='lesion', vars=session.lesion))

    #    
    #    components = [(INPUT(_type="checkbox", _name=region), region) for region in fullnames]
    #    
    #    form2 = FORM(*components) #,INPUT(_type='submit', _value='Go')) # _onclick=URL(r=request,f='compare_view',vars=vars)))
    #    # *[LI(INPUT(_type="checkbox", _name=region), region) for region in fullnames]

    #                
    #    for fullname in fullnames:
    #        if fullname in request.vars:
    #            lesion_regions.append(fullname)
    #return dict(form=form, fullnames=fullnames, vars=vars)
    return dict(form=form, form2=form2, fullnames_filepath=fullnames_filepath)
    
def regionnames():
    network = db(db.upload_data.id==request.vars.network_name).select().first()
    fullnames_filepath = os.path.join(request.folder,'uploads',network.region_names_full_file)
    fullnames = file_reader(fullnames_filepath,1)
    return fullnames

@cache(request.env.path_info, time_expire=300, cache_model=cache.ram)
def single_view():
    matrix_abs = False # Take absolute value of connectivity matrix
    matrix_minshift = True # Shift matrix values so negative weights become weakest positive weights
    network = db(db.upload_data.network_name==request.vars.network_name).select().first()
    session.network = network
    matfilename = network.connectivity_matrix_file
    matfilepath = os.path.join(request.folder,'uploads',matfilename)   
    matfilepaths = []
    matfilepaths.append(matfilepath)
    centersfilepaths = []
    centersfilepaths.append(os.path.join(request.folder,'uploads',network.region_xyz_centers_file))
    namesfilepaths = []
    namesfilepaths.append(os.path.join(request.folder,'uploads',network.region_names_full_file))
    fullnames = []
    abbrevfilepath = []
    abbrevfilepath.append(os.path.join(request.folder,'uploads',network.region_names_abbrev_file))
    weight = []
    weight.append(request.vars.weight)
    edge_density = []
    edge_density.append(request.vars.edge_density)
    names = []
    names.append(file_reader(abbrevfilepath[0],1))
    graph = []
    cmat = []
    centers = []
    names_dict = []
    num_regions = []
    orientation = request.vars.orientation
    
    fullnames_file = namesfilepaths[0]
    fullnames.append(file_reader(fullnames_file,1))

    threshold_pct = int(edge_density[0])
    if weight[0] == 'Binary':
        binarize = 1
    else:
        binarize = 0
    edge_interval_pct = 10
    connectmat_file = matfilepaths[0]
    centers_file = centersfilepaths[0]
    names_file = abbrevfilepath[0]
    cmat_pre = file_reader(connectmat_file)       
    cur_cmat = np.array(cmat_pre)
    cur_cmat[np.isinf(cur_cmat)] = 0
    cur_cmat[np.isnan(cur_cmat)] = 0
    raw_cmat = cur_cmat
    num_regions.append(len(cur_cmat))
    raw_density = (len(np.nonzero(cur_cmat)[0])/float(cur_cmat.shape[0]*(cur_cmat.shape[1]-1)))*100
    chosen_density = threshold_pct
    
    # transformation of weights
    if matrix_abs:
        cur_cmat = abs(raw_cmat)
    if matrix_minshift:
        min = np.amin(cur_cmat)
        if min >= 0:
            pass
        else:
            cur_cmat = cur_cmat - min
    if threshold_pct:
        #thresh = scipy.stats.scoreatpercentile(cur_cmat.ravel(),100-threshold_pct)
        thresh = my_scoreatpercentile(cur_cmat, threshold_pct)
        cmat_thresh = cur_cmat*(cur_cmat > thresh)
    else:
        cmat_thresh = cur_cmat
    raw_cmat_thresh = raw_cmat * (cmat_thresh > 0)
    
    edge_attributes = Storage({}) # web2py Storage class extends a dictionary
    # find mean and stddev of edge weight for all edges in raw matrix, before any potential weight transformation
    # assumes symmetric matrix, only looks at upper triangle
    edge_attributes['cur_cmat_nz_mean'] = raw_cmat[np.nonzero(raw_cmat)].mean()
    edge_attributes['cur_cmat_nz_std'] = raw_cmat[np.nonzero(raw_cmat)].std()
    
    # find mean and stddev of edge weight for edges above threshold
    edge_attributes['cmat_thresh_nz_mean'] = raw_cmat_thresh[np.nonzero(raw_cmat_thresh)].mean()
    edge_attributes['cmat_thresh_nz_std'] = raw_cmat_thresh[np.nonzero(raw_cmat_thresh)].std()
    
    if binarize:
        cmat_thresh = 1*(cmat_thresh != 0)
    cmat.append(cmat_thresh)
    G = nx.from_numpy_matrix(cmat_thresh)
    graph.append(G)
    
    if names_file:
        cur_names = file_reader(names_file,1)
        cur_names_dict={}
        for i in range(len(cur_names)):
            cur_names_dict[i] = cur_names[i]
    names_dict.append(cur_names_dict)
    
    cur_centers = core.file_reader(centers_file)
    centersa = np.array(cur_centers)

    # find mean and stddev of euclidean distance for edges above threshold
    eucdistmat = core.euclidean_distance(centersa)
    cmat_raw_mask = (cur_cmat > 0) * 1
    cmat_thresh_mask = (cmat_thresh > 0) * 1
    
    eucdistmat_masked_raw = eucdistmat * cmat_raw_mask
    eucdistmat_masked_thresh = eucdistmat * cmat_thresh_mask

    edge_attributes['eucdist_mean'] = eucdistmat_masked_thresh[np.nonzero(eucdistmat_masked_thresh)].mean()
    edge_attributes['eucdist_std'] = eucdistmat_masked_thresh[np.nonzero(eucdistmat_masked_thresh)].std()   
    
    # find mean and stddev of euclidean distance for all edges in raw matrix
    edge_attributes['eucdist_unthresh_mean'] = eucdistmat_masked_raw[np.nonzero(np.triu(eucdistmat_masked_raw, 1))].mean()
    edge_attributes['eucdist_unthresh_std'] = eucdistmat_masked_raw[np.nonzero(np.triu(eucdistmat_masked_raw, 1))].std()
    
    if orientation == 'Sagittal':
        centersxy = centersa[:,1:]
    elif orientation == 'Coronal':
        centersxy = np.column_stack((centersa[:,0],centersa[:,2]))
    else:
        centersxy = centersa[:,0:2]
    centers.append(centersxy)
    
    node_metrics={}
    if binarize:
        bc = nx.betweenness_centrality(G)
        cc = nx.clustering(G)
        ccs = []
        cur_ccs = np.array([cc[x] for x in cc])
        ccs.append(cur_ccs)
        mcc = np.mean(cur_ccs)
        try:
            cpl = nx.average_shortest_path_length(G)
        except nx.NetworkXError as e:
            cpl = 'Inf'
        sigma,lambda_cc,gamma_cpl = core.nx_small_worldness(G, True, 5000, 1, mcc, cpl)
        e_glob, cur_e_regs = core.global_efficiency(G, regional=True, weight=False)
        cur_e_regs = np.array(cur_e_regs)
        #cur_e_locs = core.local_efficiency(G, weight=False)
        cur_part_coefs = core.participation_coefficient(G)
    else:
        inv_cmat_thresh = 1./cmat_thresh
        inv_cmat_thresh[np.isinf(inv_cmat_thresh)] = 0
        G_inv = nx.from_numpy_matrix(inv_cmat_thresh)
        bc = nx.betweenness_centrality(G_inv, weight='weight')
        cc = nx.clustering(G_inv, weight='weight')
        ccs = []
        cur_ccs = np.array([cc[x] for x in cc])
        ccs.append(cur_ccs)
        mcc = np.mean(cur_ccs)
        try:
            cpl = nx.average_shortest_path_length(G_inv, weight='weight')
        except nx.NetworkXError as e:
            cpl = 'Inf'
        sigma,lambda_cc,gamma_cpl = core.nx_small_worldness(G, False, 5000, 1, mcc, cpl)
        e_glob, cur_e_regs = core.global_efficiency(G_inv, regional=True, weight=True)
        cur_e_regs = np.array(cur_e_regs)
        #cur_e_locs = core.local_efficiency(G_inv, weight=True)
        cur_part_coefs = core.participation_coefficient(G, weighted_edges=True)
    cur_part_coefs[np.isnan(cur_part_coefs)] = 0
    e_regs = []
    e_regs.append(cur_e_regs)
    #e_locs = []
    #e_locs.append(cur_e_locs)
    cur_bcs = np.array([bc[x] for x in bc])
    bcs = []
    bcs.append(cur_bcs)
    deg = nx.degree(G)
    cur_degs = np.array([deg[x] for x in deg])
    degs = []
    degs.append(cur_degs)
    num_comp = nx.number_connected_components(G)
    partition = community.best_partition(G)
    partition_list = []
    partitions = []
    for count in range(len(partition)):
        partition_list.append(partition[count]) # partition is a dictionary
    partitions.append(partition_list)
    q = community.modularity(partition, G)
    part_coefs = []
    part_coefs.append(cur_part_coefs)
    
    # Generate images
    #cmatimage = URL(r=request,f='image_mat',args=[0])
    
    cmat_filepath = os.path.join(request.folder,'static', network.network_name \
                 + '_' + weight[0] + '_' + str(chosen_density) \
                 + '_' + 'cmat' + '.png')
    degreedist_filepath = os.path.join(request.folder,'static', network.network_name \
                 + '_' + weight[0] + '_' + str(chosen_density) \
                 + '_' + 'degreedist' + '.png')
    net1_filepath = os.path.join(request.folder,'static', network.network_name \
                 + '_' + weight[0] + '_' + str(chosen_density) \
                 + '_' + orientation + '_' + 'net1' + '.png')
    bar1_filepath = os.path.join(request.folder,'static', network.network_name \
                 + '_' + weight[0] + '_' + str(chosen_density) \
                 + '_' + 'bar1' + '.png')            
    net2_filepath = os.path.join(request.folder,'static', network.network_name \
                 + '_' + weight[0] + '_' + str(chosen_density) \
                 + '_' + orientation + '_' + 'net2' + '.png')            
    bar2_filepath = os.path.join(request.folder,'static', network.network_name \
                 + '_' + weight[0] + '_' + str(chosen_density) \
                 + '_' + 'bar2' + '.png')
    net3_filepath = os.path.join(request.folder,'static', network.network_name \
                 + '_' + weight[0] + '_' + str(chosen_density) \
                 + '_' + orientation + '_' + 'net3' + '.png')
    bar3_filepath = os.path.join(request.folder,'static', network.network_name \
                 + '_' + weight[0] + '_' + str(chosen_density) \
                 + '_' + 'bar3' + '.png')
    ereg_network_filepath = os.path.join(request.folder,'static', network.network_name \
                 + '_' + weight[0] + '_' + str(chosen_density) \
                 + '_' + orientation + '_' + 'ereg_network' + '.png')    
    ereg_hist_filepath = os.path.join(request.folder,'static', network.network_name \
                 + '_' + weight[0] + '_' + str(chosen_density) \
                 + '_' + 'ereg_hist' + '.png')
    partcoef_network_filepath = os.path.join(request.folder,'static', network.network_name \
                 + '_' + weight[0] + '_' + str(chosen_density) \
                 + '_' + orientation + '_' + 'partcoef_network' + '.png')    
    partcoef_hist_filepath = os.path.join(request.folder,'static', network.network_name \
                 + '_' + weight[0] + '_' + str(chosen_density) \
                 + '_' + 'partcoef_hist' + '.png')
    module_filepath = os.path.join(request.folder,'static', network.network_name \
                 + '_' + weight[0] + '_' + str(chosen_density) \
                 + '_' + 'module' + '.png')
    spring_filepath = os.path.join(request.folder,'static', network.network_name \
                 + '_' + weight[0] + '_' + str(chosen_density) \
                 + '_' + 'spring' + '.png')
    
    filepaths = [cmat_filepath, degreedist_filepath, net1_filepath, bar1_filepath,
                 net2_filepath, bar2_filepath, net3_filepath, bar3_filepath,
                 ereg_network_filepath, ereg_hist_filepath, partcoef_network_filepath, partcoef_hist_filepath,
                 module_filepath, spring_filepath]
    
    # Make image files if they do not yet exist
    if not all([os.path.exists(file) for file in filepaths]):
        cmatimage = image_mat(0, cmat_filepath, cmat).split('/')[-1]
        degreedistimage = degree_dist(0, degreedist_filepath, graph).split('/')[-1]
        netimage1 = plot_matrix_2d(0, 's', net1_filepath, edge_density, weight, \
                                   abbrevfilepath, names_dict, graph, centers, \
                                   cmat, degs).split('/')[-1]
        barimage1 = plot_bars(0, bar1_filepath, num_regions, degs, names).split('/')[-1]
        netimage2 = plot_matrix_2d(0, 'bc', net2_filepath, edge_density, weight, \
                                   abbrevfilepath, names_dict, graph, centers, \
                                   cmat, bcs).split('/')[-1]
        barimage2 = plot_bars(0, bar2_filepath, num_regions, bcs, names).split('/')[-1]
        netimage3 = plot_matrix_2d(0, 'cc', net3_filepath, edge_density, weight, \
                                   abbrevfilepath, names_dict, graph, centers, \
                                   cmat, ccs).split('/')[-1]
        barimage3 = plot_bars(0, bar3_filepath, num_regions, ccs, names).split('/')[-1]
        ereg_hist = plot_bars(0, ereg_hist_filepath, num_regions, e_regs, names).split('/')[-1]
        ereg_network = plot_matrix_2d(0, 'ereg', ereg_network_filepath, edge_density, weight, \
                                   abbrevfilepath, names_dict, graph, centers, \
                                   cmat, e_regs).split('/')[-1]
        partcoef_hist = plot_bars(0, partcoef_hist_filepath, num_regions, part_coefs, names).split('/')[-1]
        partcoef_network = plot_matrix_2d(0, 'part_coef', partcoef_network_filepath, edge_density, weight, \
                                   abbrevfilepath, names_dict, graph, centers, \
                                   cmat, part_coefs).split('/')[-1]
        
        moduleimage = plot_modules(0, 0, module_filepath, matfilepaths, weight, names_dict, centers, graph, cmat).split('/')[-1]
        springimage = plot_modules(0, 1, spring_filepath, matfilepaths, weight, names_dict, centers, graph, cmat).split('/')[-1]
    else:
        cmatimage = cmat_filepath.split('/')[-1]
        degreedistimage = degreedist_filepath.split('/')[-1]
        netimage1 = net1_filepath.split('/')[-1]
        barimage1 = bar1_filepath.split('/')[-1]
        netimage2 = net2_filepath.split('/')[-1]
        barimage2 = bar2_filepath.split('/')[-1]
        netimage3 = net3_filepath.split('/')[-1]
        barimage3 = bar3_filepath.split('/')[-1]
        ereg_hist = ereg_hist_filepath.split('/')[-1]
        ereg_network = ereg_network_filepath.split('/')[-1]
        partcoef_hist = partcoef_hist_filepath.split('/')[-1]
        partcoef_network = partcoef_network_filepath.split('/')[-1]
        # Need to run this one every time in order to get module colors for webgl
        moduleimage = plot_modules(0, 0, module_filepath, matfilepaths, weight, names_dict, centers, graph, cmat).split('/')[-1]
        springimage = spring_filepath.split('/')[-1]
    import time
    t = time.ctime()
    
    metrics = []
    for count,region in enumerate(fullnames[0]):
        metrics.append([region,degs[0][count],
                        ccs[0][count],\
                        bcs[0][count],\
                        partitions[0][count],
                        e_regs[0][count],
                        part_coefs[0][count]]) #,e_locs[0][count]])
    session.metrics = metrics

    matfilepath = os.path.join(request.folder,'uploads',matfilename)
    centersfilepath = os.path.join(request.folder,'uploads',network.region_xyz_centers_file)
    webgl_url = plot_matrix_webgl(matfilepath,centersfilepath, network, weight=weight[0], chosen_density=chosen_density)
    
    d = dict(network=network,raw_density=raw_density,chosen_density=chosen_density,
                cmatimage=cmatimage,degreedistimage=degreedistimage,
                netimage1=netimage1,barimage1=barimage1,netimage2=netimage2,barimage2=barimage2,
                netimage3=netimage3,barimage3=barimage3,moduleimage=moduleimage,springimage=springimage,
                ereg_hist=ereg_hist,ereg_network=ereg_network,
                partcoef_hist=partcoef_hist,partcoef_network=partcoef_network,
                a=edge_density[0],b=abbrevfilepath[0],
                cpl=cpl,mcc=mcc,num_comp=num_comp,sigma=sigma,lambda_cc=lambda_cc,gamma_cpl=gamma_cpl,
                q=q, time=t, weight=weight[0], webgl_url=webgl_url, num_regions=num_regions,
                edge_attributes=edge_attributes)
    
    report_filename = pdfreport(0, weight, network, chosen_density, raw_density, cpl,\
                                mcc, q, num_comp, sigma, lambda_cc, gamma_cpl,\
                                cmatimage, degreedistimage, barimage1, netimage1,\
                                barimage2, netimage2, barimage3, netimage3,\
                                moduleimage, springimage, metrics)
    
    d["report_filename"] = report_filename
    d["metrics_filename"] = metrics_text(metrics, network, weight[0], chosen_density)
    d["e_glob"] = e_glob
    #return response.render(d)
    
    row = db(db.analysis_results.network_name==network.id).select().first()
    db(db.analysis_results.network_name==network.id).update(num_times_analyzed=row.num_times_analyzed+1)
    now = datetime.datetime.now()
    db(db.analysis_results.network_name==network.id).update(last_time_analyzed=now.strftime("%Y-%m-%d %H:%M"))

    return d

def file_reader(infile,text=0):
    fin = open(infile,'rU')
    values = []
    if text==1:
        for line in fin:
            pos = line.rstrip()
            values.append(pos)
    else:
        for line in fin:
            pos = line.rstrip().split()
            values.append(map(float, pos))
    fin.close()
    return values

def double():   
    if auth.user:
        rows_allowed = db((db.upload_data.share=='Public') | (db.upload_data.email==auth.user.email)).select()
        study_names = list(set([row.study_name for row in rows_allowed]))
        if request.vars.study_name_1 or request.vars.study_name_2:
            study_networks = [row.network_name for row in rows_allowed if row.study_name==request.vars.study_name]
        else:
            study_networks = ''
        form = SQLFORM.factory(
            Field('study_name_1', requires=IS_IN_SET(study_names),
                  default = (lambda x: x if x else '')(request.vars.study_name_cur)),
            Field('network_1_name',
                  #db.upload_data,
                  #requires=IS_IN_DB(db((db.upload_data.share=='Public') | (db.upload_data.email==auth.user.email)),
                  #                      'upload_data.id',
                  #                      '%(network_name)s')),
                  requires=IS_IN_SET(study_networks)),
            Field('study_name_2', requires=IS_IN_SET(study_names),
                  default = (lambda x: x if x else '')(request.vars.study_name_cur)),
            Field('network_2_name',
                  #db.upload_data,
                  #requires=IS_IN_DB(db((db.upload_data.share=='Public') | (db.upload_data.email==auth.user.email)),
                  #                      'upload_data.id',
                  #                      '%(network_name)s')
                  requires=IS_IN_SET(study_networks)),
            Field('weight_1',requires=IS_IN_SET([('Binary','Binary'),('Weighted','Weighted')]), default='Binary'),
            Field('weight_2',requires=IS_IN_SET([('Binary','Binary'),('Weighted','Weighted')]), default='Binary'),
            Field('edge_density_1',default=20,requires=IS_INT_IN_RANGE(0, 101, error_message=T('Must be integer in range 0-100'))),
            Field('edge_density_2',default=20,requires=IS_INT_IN_RANGE(0, 101, error_message=T('Must be integer in range 0-100'))),
            Field('orientation',requires=IS_IN_SET([('Axial','Axial'),('Sagittal','Sagittal'),('Coronal','Coronal')]), default='Axial'),
            submit_button='Analyze')
    else:
        rows_allowed = db((db.upload_data.share=='Public')).select()
        study_names = list(set([row.study_name for row in rows_allowed]))
        if request.vars.study_name:
            study_networks = [row.network_name for row in rows_allowed if row.study_name==request.vars.study_name]
        else:
            study_networks = ''
        form = SQLFORM.factory(
            Field('study_name_1', requires=IS_IN_SET(study_names),
                  default = (lambda x: x if x else '')(request.vars.study_name_cur)),
            Field('network_1_name',
                  requires=IS_IN_SET(study_networks)),
            Field('study_name_2', requires=IS_IN_SET(study_names),
                  default = (lambda x: x if x else '')(request.vars.study_name_cur)),
            Field('network_2_name',
                  requires=IS_IN_SET(study_networks)),
            Field('weight_1',requires=IS_IN_SET([('Binary','Binary'),('Weighted','Weighted')]), default='Binary'),
            Field('weight_2',requires=IS_IN_SET([('Binary','Binary'),('Weighted','Weighted')]), default='Binary'),
            Field('edge_density_1',default=20,requires=IS_INT_IN_RANGE(0, 101, error_message=T('Must be integer in range 0-100'))),
            Field('edge_density_2',default=20,requires=IS_INT_IN_RANGE(0, 101, error_message=T('Must be integer in range 0-100'))),
            Field('orientation',requires=IS_IN_SET([('Axial','Axial'),('Sagittal','Sagittal'),('Coronal','Coronal')]), default='Axial'),
            submit_button='Analyze')
    form[0][0].insert(-1,"Choose from the list of available studies")
    form[0][1].insert(-1,"Choose from the list of shared networks for the first chosen study")
    form[0][2].insert(-1,"Choose from the list of available studies")
    form[0][3].insert(-1,"Choose from the list of shared networks for the second chosen study")
    form[0][4].insert(-1,"Choose a weighting scheme for the network, either weighted or binarized edges")
    form[0][5].insert(-1,"Choose a weighting scheme for the network, either weighted or binarized edges")
    form[0][6].insert(-1,"Choose the percent of connections to keep (0-100%)")
    form[0][7].insert(-1,"Choose the percent of connections to keep (0-100%)")
    form[0][8].insert(-1,"Choose the orientation from which to view the network")
    if form.accepts(request.vars, session):
        cache.ram.clear()
        redirect(URL(r=request,f='compare_view',vars=form.vars))
    return dict(form=form)
    
def about():
    return dict()

def image_mat(num, path, cmat):
    plt.figure(figsize=(5,5))
    plt.imshow(cmat[num],interpolation='nearest')
    plt.colorbar()
    savefig(path)
    plt.close()
    return path
    
def degree_dist(num, path, graph):
    deg = nx.degree(graph[num])
    degs = np.array([deg[x] for x in deg])
    plt.figure(figsize=(5,5))
    plt.hist(degs,10)
    plt.title("Node Degree Histogram")
    plt.ylabel("Degree")
    plt.xlabel("Number")
    savefig(path)
    plt.close()
    return path
    
#def plot_bars(num, node_metric, path, num_regions, bcs, degs, ccs, names):
def plot_bars(num, path, num_regions, metric, names):
    """
    Plot histogram of values of given regional metric in descending order
    """
    num_regions = num_regions[num]
    matplotlib.rcParams['font.size'] = 6
    val = metric[num]
    #if node_metric == 's':
    #    val = degs[num]
    #elif node_metric == 'bc':
    #    val = bcs[num]
    #elif node_metric == 'cc':
    #    val = ccs[num]
    val_sort_indices=np.argsort(val)
    val_sort=[val[x] for x in val_sort_indices]
    val_names_sort=[names[num][x] for x in val_sort_indices]  
    pos = np.arange(num_regions)*2
    figheight = num_regions/12.
    plt.figure(figsize = (5, figheight))    
    plt.barh(pos, val_sort)
    plt.yticks(pos, val_names_sort)
    plt.autoscale(tight=True)
    #plt.show()
    savefig(path)
    plt.close()
    matplotlib.rcParams['font.size'] = 12
    return path
    
def plot_modules(num, spring, path, matfilepaths, weight, names_dict, centers, graph, cmat):
    import matplotlib.cm
    import matplotlib.colors as colors
    
    connectmat_file = matfilepaths[num]

    names_dict = names_dict[num]    
    centersxy = centers[num]
    G = graph[num]
    cmat_thresh = cmat[num]
    
    if weight[num] == 'Binary':
        binarize = 1
        weight_edges = False
    else:
        binarize = 0
        weight_edges = True
    alpha = .5
    edge_interval_pct = 10
    # orientation; get proper aspect ratio first
    xy_size = (5,5)
    if not spring:
        basesize = 5
        xr = abs(min(centersxy[:,0])) + abs(max(centersxy[:,0]))
        yr = abs(min(centersxy[:,1])) + abs(max(centersxy[:,1]))
    
        xy_ratio = max(xr,yr)/min(xr,yr)
        if xr > yr:
            xy_size = (basesize*xy_ratio, basesize)
        else:
            xy_size = (basesize, basesize*xy_ratio)

    plt.figure(figsize=xy_size)
    
    partition = community.best_partition(G)

    size = float(len(set(partition.values())))
    pos = nx.spring_layout(G)
    N=np.array(range(int(size)))
    
    fracs = N.astype(float)/N.max()
    norm = colors.normalize(fracs.min(), fracs.max())
    
    module_colors = [0]*len(names_dict)
    count = 0.
    for com in set(partition.values()) :
        count = count + 1.
        list_nodes = [nodes for nodes in partition.keys()
                                    if partition[nodes] == com]
        if spring == 1:
            nx.draw_networkx_nodes(G, pos, list_nodes, node_size = 200,
                                    node_color = matplotlib.cm.jet(norm(fracs[count-1])))
        else:
            nx.draw_networkx_nodes(G, centersxy, list_nodes, node_size = 200,
                                    node_color = matplotlib.cm.jet(norm(fracs[count-1])))
        rgb = matplotlib.cm.jet(norm(fracs[count-1]))[0:3]
        rgb_255 = tuple([int(a*255) for a in rgb])
        hex = '#%02x%02x%02x' % tuple(rgb_255)
        for node in list_nodes:
            module_colors[node] = hex
    session.module_colors = module_colors

    if weight_edges:
        edges = []
        nonzero_edges = cmat_thresh[np.nonzero(cmat_thresh)] # all nonzero edges
        percentiles = [my_scoreatpercentile(nonzero_edges, 100-x) for x in range(0,101,edge_interval_pct)]
        for i in range(len(percentiles)-1):
            alpha_val = .1 + (i / 20.0) # edges in first percentile have alpha=0
            thresh_low = percentiles[i]
            thresh_high = percentiles[i+1]
            edges.append([(u,v) for (u,v,d) in G.edges(data=True) if thresh_low < d['weight'] < thresh_high])
            if spring == 1:
                nx.draw_networkx_edges(G,pos,edgelist=edges[i],width=i/1.9,alpha=alpha_val,edge_color='k')
            else:
                nx.draw_networkx_edges(G,centersxy,edgelist=edges[i],width=i/1.9,alpha=alpha_val,edge_color='k')
    else:
        if spring == 1:
            nx.draw_networkx_edges(G, pos, width=1, alpha=alpha)
        else:
            nx.draw_networkx_edges(G, centersxy, width=1, alpha=alpha)
    
    if spring == 1:
        nx.draw_networkx_labels(G, pos, labels=names_dict, font_size=10)
    else:
        nx.draw_networkx_labels(G, centersxy, labels=names_dict, font_size=10)
    plt.autoscale(tight=True)
    savefig(path)
    plt.close()
    return path

def plot_matrix_2d(num, node_metric, path, edge_density, weight, abbrevfilepath, \
                   names_dict, graph, centers, cmat, metric):
    num = int(num)
    grp_metrics = 0
    threshold_pct = int(edge_density[num])
    if weight[num] == 'Binary':
        binarize = 1
        weight_edges = False
    else:
        binarize = 0
        weight_edges = True
    alpha = .5
    edge_interval_pct = 10
    names_file = abbrevfilepath[num]
    names_dict = names_dict[num]
    G = graph[num]
    centersxy = centers[num]
    cmat_thresh = cmat[num]
    node_colors = {}     # define colors for each metric
    node_colors['s'] = 'orange'
    node_colors['cc'] = 'aqua'
    node_colors['bc'] = 'chartreuse'
    node_colors['ereg'] = 'magenta'
    node_colors['part_coef'] = 'yellow'
    node_metrics={}
    if grp_metrics: # regional metrics caclulated elsewhere, loaded in
        metrics = np.array(file_reader(grp_metrics))
        cols = np.shape(metrics)[1]
        for i in range(cols):
            colmean = np.mean(metrics[:,i])
            colscale = 300 / colmean
            metrics[:,i] = metrics[:,i] * colscale # make node mean radius 300
        node_metrics['s'] = metrics[:,0] # strength
        node_metrics['cc'] = metrics[:,1] # clustering coefficient
        node_metrics['bc'] = metrics[:,2] # betweenness centrality
        node_metrics['eloc'] = metrics[:,3] # local efficiency
        node_metrics['ereg'] = metrics[:,4] # regional efficiency
    else: # otherwise, calculate them locally (networkx measures somewhat incomplete)
        #cur_bcs = bcs[num]
        #bcscale = 200 / np.mean(cur_bcs)
        #cur_bcs = cur_bcs * bcscale
        #node_metrics['bc'] = cur_bcs
        #cur_degs = degs[num]
        #degscale = 200 / np.mean(cur_degs)
        #cur_degs = cur_degs * degscale
        #node_metrics['s'] = cur_degs
        #cur_ccs = ccs[num]
        #ccscale = 200 / np.mean(cur_ccs)
        #cur_ccs = cur_ccs * ccscale
        #node_metrics['cc'] = cur_ccs
        cur_metric = metric[num]
        metric_scaled = 200 / np.mean(cur_metric)
        cur_metric = cur_metric * metric_scaled
    
    # orientation; get proper aspect ratio first
    basesize = 5
    xr = abs(min(centersxy[:,0])) + abs(max(centersxy[:,0]))
    yr = abs(min(centersxy[:,1])) + abs(max(centersxy[:,1]))

    xy_ratio = max(xr,yr)/min(xr,yr)
    if xr > yr:
        xy_size = (basesize*xy_ratio, basesize)
    else:
        xy_size = (basesize, basesize*xy_ratio)

    plt.figure(figsize=xy_size)
    
    #print cur_metric
    #print node_colors[node_metric]
    nx.draw_networkx_nodes(G, centersxy, node_size=cur_metric, node_color=node_colors[node_metric])
    if names_file:
        nx.draw_networkx_labels(G,centersxy,labels=names_dict,font_size=10)
    if weight_edges:
        edges = []
        #percentiles = [scipy.stats.scoreatpercentile(ma.ravel(),x) for x in range(0,101,edge_interval_pct)]
        nonzero_edges = cmat_thresh[np.nonzero(cmat_thresh)] # all nonzero edges
        #percentiles = np.histogram(nonzero_edges,10)[1] # equally binned histogram of non-zero edges
        #percentiles = [scipy.stats.scoreatpercentile(nonzero_edges,x) for x in range(0,101,edge_interval_pct)]
        percentiles = [my_scoreatpercentile(nonzero_edges, 100-x) for x in range(0,101,edge_interval_pct)]
        for i in range(len(percentiles)-1):
            alpha_val = .1 + (i / 20.0) # edges in first percentile have alpha=0
            thresh_low = percentiles[i]
            thresh_high = percentiles[i+1]
            edges.append([(u,v) for (u,v,d) in G.edges(data=True) if thresh_low < d['weight'] < thresh_high])
            nx.draw_networkx_edges(G,centersxy,edgelist=edges[i],width=i/1.9,alpha=alpha_val,edge_color='k')
        plt.autoscale(tight=True)
        #plt.show()
        savefig(path)
        plt.close()
        return path
    else:
        nx.draw_networkx_edges(G,centersxy,width=1,alpha=alpha,edge_color='k')
        plt.autoscale(tight=True)
        #plt.show()
        savefig(path)
        plt.close()
        return path

@cache(request.env.path_info, time_expire=300, cache_model=cache.ram)
def compare_view():
    matrix_abs = False # Take absolute value of connectivity matrix
    matrix_minshift = True # Shift matrix values so negative weights become weakest positive weights
    network1 = db(db.upload_data.network_name==request.vars.network_1_name).select().first()
    session.network1 = network1
    network2 = db(db.upload_data.network_name==request.vars.network_2_name).select().first()
    session.network2 = network2
    matfilename1 = network1.connectivity_matrix_file
    matfilename2 = network2.connectivity_matrix_file
    matfilepath1 = os.path.join(request.folder,'uploads',matfilename1)
    matfilepath2 = os.path.join(request.folder,'uploads',matfilename2)
    
    matfilepaths = []
    matfilepaths.append(matfilepath1)
    matfilepaths.append(matfilepath2)
    centersfilepaths = []
    centersfilepaths.append(os.path.join(request.folder,'uploads',network1.region_xyz_centers_file))
    centersfilepaths.append(os.path.join(request.folder,'uploads',network2.region_xyz_centers_file))
    namesfilepaths = []
    namesfilepaths.append(os.path.join(request.folder,'uploads',network1.region_names_full_file))
    namesfilepaths.append(os.path.join(request.folder,'uploads',network2.region_names_full_file))
    abbrevfilepaths = []
    abbrevfilepaths.append(os.path.join(request.folder,'uploads',network1.region_names_abbrev_file))
    abbrevfilepaths.append(os.path.join(request.folder,'uploads',network2.region_names_abbrev_file))
    weights = []
    weights.append(request.vars.weight_1)
    weights.append(request.vars.weight_2)
    edge_densities = []
    edge_densities.append(request.vars.edge_density_1)
    edge_densities.append(request.vars.edge_density_2)
    names = []
    names.append(file_reader(abbrevfilepaths[0],1))
    names.append(file_reader(abbrevfilepaths[1],1))
    graphs = []
    cmats = []
    centers = []
    names_dicts = []
    num_regions = []
    bcs = []
    degs = []
    ccs = []
    fullnames = []
    partitions = []
    orientation = request.vars.orientation
    e_regs = []
    #e_locs = []
    part_coefs = []
    
    raw_densities = []
    chosen_densities = []
    cpl = []
    mcc = []
    num_comp = []
    sigma = []
    lambda_cc = []
    gamma_cpl = []
    q = []
    e_glob = []
    edge_attributes_list = []
    
    for compare_num in range(2):
        fullnames_file = namesfilepaths[compare_num]
        cur_fullnames = file_reader(fullnames_file,1)
        fullnames.append(cur_fullnames)
        
        threshold_pct = int(edge_densities[compare_num])
        if weights[compare_num] == 'Binary':
            binarize = 1
        else:
            binarize = 0
        edge_interval_pct=10
        connectmat_file = matfilepaths[compare_num]
        centers_file = centersfilepaths[compare_num]
        names_file = abbrevfilepaths[compare_num]
        cmat_pre = file_reader(connectmat_file)
        cur_cmat = np.array(cmat_pre)
        cur_cmat[np.isinf(cur_cmat)] = 0
        cur_cmat[np.isnan(cur_cmat)] = 0
        raw_cmat = cur_cmat
        num_regions.append(len(cur_cmat))
        raw_densities.append((len(np.nonzero(cur_cmat)[0])/float(cur_cmat.shape[0]*(cur_cmat.shape[1]-1)))*100)
        chosen_densities.append(threshold_pct)
        # set lesioned nodes to 0
        if compare_num == 1 and request.args:
            for r in session.lesion['region_nums']:
            #for r in eval(request.vars["region_nums"]):
                cur_cmat[r,:] = 0
                cur_cmat[:,r] = 0
        if matrix_abs:
            cur_cmat = abs(cur_cmat)
        if matrix_minshift:
            min = np.amin(cur_cmat)
            if min >= 0:
                pass
            else:
                cur_cmat = cur_cmat - min
        if threshold_pct:
            #thresh = scipy.stats.scoreatpercentile(cur_cmat.ravel(),100-threshold_pct)
            thresh = my_scoreatpercentile(cur_cmat, threshold_pct)
            cmat_thresh = cur_cmat*(cur_cmat > thresh)
            cur_cmat_thresh = cur_cmat*(cur_cmat > thresh) # unbinarized
        else:
            cur_cmat_thresh = cur_cmat
        if binarize:
            cmat_thresh = 1*(cmat_thresh != 0) # binarized
        cmats.append(cmat_thresh)
        G = nx.from_numpy_matrix(cmat_thresh)
        graphs.append(G)

        raw_cmat_thresh = raw_cmat * (cur_cmat_thresh > 0)
        
        edge_attributes = Storage({}) # web2py Storage class extends a dictionary
        edge_attributes_list.append(edge_attributes)
        # find mean and stddev of edge weight for all edges in raw matrix, before any potential weight transformation
        # assumes symmetric matrix, only looks at upper triangle
        edge_attributes['cur_cmat_nz_mean'] = raw_cmat[np.nonzero(raw_cmat)].mean()
        edge_attributes['cur_cmat_nz_std'] = raw_cmat[np.nonzero(raw_cmat)].std()
        
        # find mean and stddev of edge weight for edges above threshold
        edge_attributes['cmat_thresh_nz_mean'] = raw_cmat_thresh[np.nonzero(raw_cmat_thresh)].mean()
        edge_attributes['cmat_thresh_nz_std'] = raw_cmat_thresh[np.nonzero(raw_cmat_thresh)].std()
    
        if names_file:
            cur_names = file_reader(names_file,1)
            cur_names_dict={}
            for i in range(len(cur_names)):
                cur_names_dict[i] = cur_names[i]
        names_dicts.append(cur_names_dict)
    
        cur_centers = core.file_reader(centers_file)
        centersa = np.array(cur_centers)
        
        if orientation == 'Sagittal':
            centersxy = centersa[:,1:]
        elif orientation == 'Coronal':
            centersxy = np.column_stack((centersa[:,0],centersa[:,2]))
        else:
            centersxy = centersa[:,0:2]
        centers.append(centersxy)
        
        # find mean and stddev of euclidean distance for edges above threshold
        eucdistmat = core.euclidean_distance(centersa)
        cmat_raw_mask = (cur_cmat > 0) * 1
        cmat_thresh_mask = (cur_cmat_thresh > 0) * 1
        
        eucdistmat_masked_raw = eucdistmat * cmat_raw_mask
        eucdistmat_masked_thresh = eucdistmat * cmat_thresh_mask
    
        edge_attributes['eucdist_mean'] = eucdistmat_masked_thresh[np.nonzero(eucdistmat_masked_thresh)].mean()
        edge_attributes['eucdist_std'] = eucdistmat_masked_thresh[np.nonzero(eucdistmat_masked_thresh)].std()   
        
        # find mean and stddev of euclidean distance for all edges in raw matrix
        edge_attributes['eucdist_unthresh_mean'] = eucdistmat_masked_raw[np.nonzero(np.triu(eucdistmat_masked_raw, 1))].mean()
        edge_attributes['eucdist_unthresh_std'] = eucdistmat_masked_raw[np.nonzero(np.triu(eucdistmat_masked_raw, 1))].std()
        
        if binarize:
            bc = nx.betweenness_centrality(G)
            cc = nx.clustering(G)
            cur_ccs = np.array([cc[x] for x in cc])
            ccs.append(cur_ccs)
            mcc.append(np.mean(cur_ccs))
            try:
                cpl.append(nx.average_shortest_path_length(G))
            except nx.NetworkXError as e:
                cpl.append('Inf')
            a, b, c = core.nx_small_worldness(G, True, 5000, 1, mcc[-1], cpl[-1])
            sigma.append(a)
            lambda_cc.append(b)
            gamma_cpl.append(c)
            cur_e_glob, cur_e_regs = core.global_efficiency(G, regional=True, weight=False)
            cur_e_regs = np.array(cur_e_regs)
            #cur_e_locs = core.local_efficiency(G, weight=False)
            cur_part_coefs = core.participation_coefficient(G)
        else:
            inv_cmat_thresh = 1 /cmat_thresh
            inv_cmat_thresh[np.isinf(inv_cmat_thresh)] = 0
            G_inv = nx.from_numpy_matrix(inv_cmat_thresh)
            bc = nx.betweenness_centrality(G_inv, weight='weight')
            cc = nx.clustering(G, weight='weight')
            cur_ccs = np.array([cc[x] for x in cc])
            ccs.append(cur_ccs)
            mcc.append(np.mean(cur_ccs))
            try:
                cpl.append(nx.average_shortest_path_length(G_inv, weight='weight'))
            except nx.NetworkXError as e:
                cpl.append('Inf')
            a, b, c = core.nx_small_worldness(G, False, 5000, 1, mcc[-1], cpl[-1])
            sigma.append(a)
            lambda_cc.append(b)
            gamma_cpl.append(c)
            cur_e_glob, cur_e_regs = core.global_efficiency(G_inv, regional=True, weight=True)
            cur_e_regs = np.array(cur_e_regs)
            #cur_e_locs = core.local_efficiency(G_inv, weight=True)
            cur_part_coefs = core.participation_coefficient(G, weighted_edges=True)
        cur_part_coefs[np.isnan(cur_part_coefs)] = 0

        #if binarize:
        #    bc = nx.betweenness_centrality(G)
        #    cc = nx.clustering(G)
        #    try:
        #        cpl.append(nx.average_shortest_path_length(G))
        #    except nx.NetworkXError as e:
        #        cpl.append('Inf')
        #    a,b,c = core.nx_small_worldness(G, True, 500, 10)
        #    sigma.append(a)
        #    lambda_cc.append(b)
        #    gamma_cpl.append(c)
        #else:
        #    bc = nx.betweenness_centrality(G,weight='weight')
        #    cc = nx.clustering(G, weighted=True)
        #    try:
        #        cpl.append(nx.average_shortest_path_length(G,weight='weight'))
        #    except nx.NetworkXError as e:
        #        cpl.append('Inf')
        #    a,b,c = core.nx_small_worldness(G, False, 500, 10)
        #    sigma.append(a)
        #    lambda_cc.append(b)
        #    gamma_cpl.append(c)
        #cur_e_glob, cur_e_regs = core.global_efficiency(G, regional=True)
        #cur_e_locs = core.local_efficiency(G)
        
        e_glob.append(cur_e_glob)
        e_regs.append(cur_e_regs)
        
        #e_locs.append(cur_e_locs)
        
        cur_bcs = np.array([bc[x] for x in bc])
        bcs.append(cur_bcs)
        deg = nx.degree_centrality(G)
        cur_degs = np.array([deg[x] for x in deg])
        degs.append(cur_degs)
        num_comp.append(nx.number_connected_components(G))
        partition = community.best_partition(G)
        partition_list = []
        for count in range(len(partition)):
            partition_list.append(partition[count]) # partition is a dictionary
        partitions.append(partition_list)
        q.append(community.modularity(partition, G))
        part_coefs.append(cur_part_coefs)

    # Generate images       
    cmat_filepath1 = os.path.join(request.folder,'static', network1.network_name \
                 + '_' + weights[0] + '_' + str(chosen_densities[0]) \
                 + '_' + 'cmat' + '.png')
    degreedist_filepath1 = os.path.join(request.folder,'static', network1.network_name \
                 + '_' + weights[0] + '_' + str(chosen_densities[0])\
                 + '_' + 'degreedist' + '.png')
    net_filepath1 = os.path.join(request.folder,'static', network1.network_name \
                 + '_' + weights[0] + '_' + str(chosen_densities[0])  + '_' + orientation \
                 + '_' + 'net1' + '.png')
    bar_filepath1 = os.path.join(request.folder,'static', network1.network_name \
                 + '_' + weights[0] + '_' + str(chosen_densities[0]) \
                 + '_' + 'bar1' + '.png')            
    net_filepath2 = os.path.join(request.folder,'static', network1.network_name \
                 + '_' + weights[0] + '_' + str(chosen_densities[0])  + '_' + orientation \
                 + '_' + 'net2' + '.png')            
    bar_filepath2 = os.path.join(request.folder,'static', network1.network_name \
                 + '_' + weights[0] + '_' + str(chosen_densities[0]) \
                 + '_' + 'bar2' + '.png')
    net_filepath3 = os.path.join(request.folder,'static', network1.network_name \
                 + '_' + weights[0] + '_' + str(chosen_densities[0]) + '_' + orientation \
                 + '_' + 'net3' + '.png')
    bar_filepath3 = os.path.join(request.folder,'static', network1.network_name \
                 + '_' + weights[0] + '_' + str(chosen_densities[0]) \
                 + '_' + 'bar3' + '.png')
    ereg_network_filepath1 = os.path.join(request.folder,'static', network1.network_name \
                 + '_' + weights[0] + '_' + str(chosen_densities[0]) \
                 + '_' + orientation + '_' + 'ereg_network' + '.png')    
    ereg_hist_filepath1 = os.path.join(request.folder,'static', network1.network_name \
                 + '_' + weights[0] + '_' + str(chosen_densities[0]) \
                 + '_' + 'ereg_hist' + '.png')
    partcoef_network_filepath1 = os.path.join(request.folder,'static', network1.network_name \
                 + '_' + weights[0] + '_' + str(chosen_densities[0]) \
                 + '_' + orientation + '_' + 'partcoef_network' + '.png')    
    partcoef_hist_filepath1 = os.path.join(request.folder,'static', network1.network_name \
                 + '_' + weights[0] + '_' + str(chosen_densities[0]) \
                 + '_' + 'partcoef_hist' + '.png')
    module_filepath1 = os.path.join(request.folder,'static', network1.network_name \
                 + '_' + weights[0] + '_' + str(chosen_densities[0])  + '_' + orientation \
                 + '_' + 'module' + '.png')
    spring_filepath1 = os.path.join(request.folder,'static', network1.network_name \
                 + '_' + weights[0] + '_' + str(chosen_densities[0]) \
                 + '_' + 'spring' + '.png')
    
    filepaths1 = [cmat_filepath1, degreedist_filepath1, net_filepath1, bar_filepath1, net_filepath2, bar_filepath2, net_filepath3, bar_filepath3, module_filepath1, spring_filepath1]
    
    if compare_num == 1 and request.args:
        cmat_filepath2 = os.path.join(request.folder,'static', network2.network_name \
                     + '_' + weights[1] + '_' + str(chosen_densities[1]) \
                     + '_' + 'lesion' + '_' + 'cmat' + '.png')
        degreedist_filepath2 = os.path.join(request.folder,'static', network2.network_name \
                     + '_' + weights[1] + '_' + str(chosen_densities[1]) \
                     + '_' + 'lesion' + '_' + 'degreedist' + '.png')
        net_filepath4 = os.path.join(request.folder,'static', network2.network_name \
                     + '_' + weights[1] + '_' + str(chosen_densities[1]) + '_' + orientation\
                     + '_' + 'lesion' + '_' + 'net1' + '.png')
        bar_filepath4 = os.path.join(request.folder,'static', network2.network_name \
                     + '_' + weights[1] + '_' + str(chosen_densities[1]) \
                     + '_' + 'lesion' + '_' + 'bar1' + '.png')            
        net_filepath5 = os.path.join(request.folder,'static', network2.network_name \
                     + '_' + weights[1] + '_' + str(chosen_densities[1]) + '_' + orientation\
                     + '_' + 'lesion' + '_' + 'net2' + '.png')            
        bar_filepath5 = os.path.join(request.folder,'static', network2.network_name \
                     + '_' + weights[1] + '_' + str(chosen_densities[1]) \
                     + '_' + 'lesion' + '_' + 'bar2' + '.png')
        net_filepath6 = os.path.join(request.folder,'static', network2.network_name \
                     + '_' + weights[1] + '_' + str(chosen_densities[1]) + '_' + orientation\
                     + '_' + 'lesion' + '_' + 'net3' + '.png')
        bar_filepath6 = os.path.join(request.folder,'static', network2.network_name \
                     + '_' + weights[1] + '_' + str(chosen_densities[1]) \
                     + '_' + 'lesion' + '_' + 'bar3' + '.png')
        ereg_network_filepath2 = os.path.join(request.folder,'static', network2.network_name \
                     + '_' + weights[1] + '_' + str(chosen_densities[1]) \
                     + '_' + orientation + '_' + 'ereg_network' + '.png')    
        ereg_hist_filepath2 = os.path.join(request.folder,'static', network2.network_name \
                     + '_' + weights[1] + '_' + str(chosen_densities[1]) \
                     + '_' + 'ereg_hist' + '.png')
        partcoef_network_filepath2 = os.path.join(request.folder,'static', network2.network_name \
                     + '_' + weights[1] + '_' + str(chosen_densities[1]) \
                     + '_' + orientation + '_' + 'partcoef_network' + '.png')    
        partcoef_hist_filepath2 = os.path.join(request.folder,'static', network2.network_name \
                     + '_' + weights[1] + '_' + str(chosen_densities[1]) \
                     + '_' + 'partcoef_hist' + '.png')
        module_filepath2 = os.path.join(request.folder,'static', network2.network_name \
                     + '_' + weights[1] + '_' + str(chosen_densities[1]) + '_' + orientation\
                     + '_' + 'lesion' + '_' + 'module' + '.png')
        spring_filepath2 = os.path.join(request.folder,'static', network2.network_name \
                     + '_' + weights[1] + '_' + str(chosen_densities[1]) \
                     + '_' + 'lesion' + '_' + 'spring' + '.png')
    else:
        cmat_filepath2 = os.path.join(request.folder,'static', network2.network_name \
                     + '_' + weights[1] + '_' + str(chosen_densities[1]) \
                     + '_' + 'cmat' + '.png')
        degreedist_filepath2 = os.path.join(request.folder,'static', network2.network_name \
                     + '_' + weights[1] + '_' + str(chosen_densities[1]) \
                     + '_' + 'degreedist' + '.png')
        net_filepath4 = os.path.join(request.folder,'static', network2.network_name \
                     + '_' + weights[1] + '_' + str(chosen_densities[1]) + '_' + orientation\
                     + '_' + 'net1' + '.png')
        bar_filepath4 = os.path.join(request.folder,'static', network2.network_name \
                     + '_' + weights[1] + '_' + str(chosen_densities[1]) \
                     + '_' + 'bar1' + '.png')            
        net_filepath5 = os.path.join(request.folder,'static', network2.network_name \
                     + '_' + weights[1] + '_' + str(chosen_densities[1]) + '_' + orientation\
                     + '_' + 'net2' + '.png')            
        bar_filepath5 = os.path.join(request.folder,'static', network2.network_name \
                     + '_' + weights[1] + '_' + str(chosen_densities[1]) \
                     + '_' + 'bar2' + '.png')
        net_filepath6 = os.path.join(request.folder,'static', network2.network_name \
                     + '_' + weights[1] + '_' + str(chosen_densities[1]) + '_' + orientation\
                     + '_' + 'net3' + '.png')
        bar_filepath6 = os.path.join(request.folder,'static', network2.network_name \
                     + '_' + weights[1] + '_' + str(chosen_densities[1]) \
                     + '_' + 'bar3' + '.png')
        ereg_network_filepath2 = os.path.join(request.folder,'static', network2.network_name \
                     + '_' + weights[1] + '_' + str(chosen_densities[1]) \
                     + '_' + orientation + '_' + 'ereg_network' + '.png')    
        ereg_hist_filepath2 = os.path.join(request.folder,'static', network2.network_name \
                     + '_' + weights[1] + '_' + str(chosen_densities[1]) \
                     + '_' + 'ereg_hist' + '.png')
        partcoef_network_filepath2 = os.path.join(request.folder,'static', network2.network_name \
                     + '_' + weights[1] + '_' + str(chosen_densities[1]) \
                     + '_' + orientation + '_' + 'partcoef_network' + '.png')    
        partcoef_hist_filepath2 = os.path.join(request.folder,'static', network2.network_name \
                     + '_' + weights[1] + '_' + str(chosen_densities[1]) \
                     + '_' + 'partcoef_hist' + '.png')
        module_filepath2 = os.path.join(request.folder,'static', network2.network_name \
                     + '_' + weights[1] + '_' + str(chosen_densities[1]) + '_' + orientation\
                     + '_' + 'module' + '.png')
        spring_filepath2 = os.path.join(request.folder,'static', network2.network_name \
                     + '_' + weights[1] + '_' + str(chosen_densities[1]) \
                     + '_' + 'spring' + '.png')
        
    filepaths2 = [cmat_filepath2, degreedist_filepath2, net_filepath4, bar_filepath4, net_filepath5, bar_filepath5, net_filepath6, bar_filepath6, module_filepath2, spring_filepath2]
    
    # Make image files if they do not yet exist
    if not all([os.path.exists(file) for file in filepaths1]):
        cmatimage1 = image_mat(0, cmat_filepath1, cmats).split('/')[-1]
        degreedistimage1 = degree_dist(0, degreedist_filepath1, graphs).split('/')[-1]
        netimage1 = plot_matrix_2d(0, 's', net_filepath1, edge_densities, weights, \
                                   abbrevfilepaths, names_dicts, graphs, centers, \
                                   cmats, degs).split('/')[-1]
        barimage1 = plot_bars(0, bar_filepath1, num_regions, degs, names).split('/')[-1]
        netimage2 = plot_matrix_2d(0, 'bc', net_filepath2, edge_densities, weights, \
                                   abbrevfilepaths, names_dicts, graphs, centers, \
                                   cmats, bcs).split('/')[-1]
        barimage2 = plot_bars(0, bar_filepath2, num_regions, bcs, names).split('/')[-1]
        netimage3 = plot_matrix_2d(0, 'cc', net_filepath3, edge_densities, weights, \
                                   abbrevfilepaths, names_dicts, graphs, centers, \
                                   cmats, ccs).split('/')[-1]
        barimage3 = plot_bars(0, bar_filepath3, num_regions, ccs, names).split('/')[-1]
        ereg_hist1 = plot_bars(0, ereg_hist_filepath1, num_regions, e_regs, names).split('/')[-1]
        ereg_network1 = plot_matrix_2d(0, 'ereg', ereg_network_filepath1, edge_densities, weights, \
                                   abbrevfilepaths, names_dicts, graphs, centers, \
                                   cmats, e_regs).split('/')[-1]
        partcoef_hist1 = plot_bars(0, partcoef_hist_filepath1, num_regions, part_coefs, names).split('/')[-1]
        partcoef_network1 = plot_matrix_2d(0, 'part_coef', partcoef_network_filepath1, edge_densities, weights, \
                                   abbrevfilepaths, names_dicts, graphs, centers, \
                                   cmats, part_coefs).split('/')[-1]
        moduleimage1 = plot_modules(0, 0, module_filepath1, matfilepaths, weights, names_dicts, centers, graphs, cmats).split('/')[-1]
        springimage1 = plot_modules(0, 1, spring_filepath1, matfilepaths, weights, names_dicts, centers, graphs, cmats).split('/')[-1]
    else:
        cmatimage1 = cmat_filepath1.split('/')[-1]
        degreedistimage1 = degreedist_filepath1.split('/')[-1]
        netimage1 = net_filepath1.split('/')[-1]
        barimage1 = bar_filepath1.split('/')[-1]
        netimage2 = net_filepath2.split('/')[-1]
        barimage2 = bar_filepath2.split('/')[-1]
        netimage3 = net_filepath3.split('/')[-1]
        barimage3 = bar_filepath3.split('/')[-1]
        ereg_hist1 = ereg_hist_filepath1.split('/')[-1]
        ereg_network1 = ereg_network_filepath1.split('/')[-1]
        partcoef_hist1 = partcoef_hist_filepath1.split('/')[-1]
        partcoef_network1 = partcoef_network_filepath1.split('/')[-1]
        moduleimage1 = module_filepath1.split('/')[-1]
        springimage1 = spring_filepath1.split('/')[-1]
    if not all([os.path.exists(file) for file in filepaths2]):
        cmatimage2 = image_mat(1, cmat_filepath2, cmats).split('/')[-1]
        degreedistimage2 = degree_dist(1, degreedist_filepath2, graphs).split('/')[-1]
        netimage4 = plot_matrix_2d(1, 's', net_filepath4, edge_densities, weights, \
                                   abbrevfilepaths, names_dicts, graphs, centers, \
                                   cmats, degs).split('/')[-1]
        barimage4 = plot_bars(1, bar_filepath4, num_regions, degs, names).split('/')[-1]
        netimage5 = plot_matrix_2d(1, 'bc', net_filepath5, edge_densities, weights, \
                                   abbrevfilepaths, names_dicts, graphs, centers, \
                                   cmats, bcs).split('/')[-1]
        barimage5 = plot_bars(1, bar_filepath5, num_regions, bcs, names).split('/')[-1]
        netimage6 = plot_matrix_2d(1, 'cc', net_filepath6, edge_densities, weights, \
                                   abbrevfilepaths, names_dicts, graphs, centers, \
                                   cmats, ccs).split('/')[-1]
        barimage6 = plot_bars(1, bar_filepath6, num_regions, ccs, names).split('/')[-1]
        ereg_hist2 = plot_bars(1, ereg_hist_filepath2, num_regions, e_regs, names).split('/')[-1]
        ereg_network2 = plot_matrix_2d(1, 'ereg', ereg_network_filepath2, edge_densities, weights, \
                                   abbrevfilepaths, names_dicts, graphs, centers, \
                                   cmats, e_regs).split('/')[-1]
        partcoef_hist2 = plot_bars(1, partcoef_hist_filepath2, num_regions, part_coefs, names).split('/')[-1]
        partcoef_network2 = plot_matrix_2d(1, 'part_coef', partcoef_network_filepath2, edge_densities, weights, \
                                   abbrevfilepaths, names_dicts, graphs, centers, \
                                   cmats, part_coefs).split('/')[-1]
        moduleimage2 = plot_modules(1, 0, module_filepath2, matfilepaths, weights, names_dicts, centers, graphs, cmats).split('/')[-1]
        springimage2 = plot_modules(1, 1, spring_filepath2, matfilepaths, weights, names_dicts, centers, graphs, cmats).split('/')[-1]
    else:
        cmatimage2 = cmat_filepath2.split('/')[-1]
        degreedistimage2 = degreedist_filepath2.split('/')[-1]
        netimage4 = net_filepath4.split('/')[-1]
        barimage4 = bar_filepath4.split('/')[-1]
        netimage5 = net_filepath5.split('/')[-1]
        barimage5 = bar_filepath5.split('/')[-1]
        netimage6 = net_filepath6.split('/')[-1]
        barimage6 = bar_filepath6.split('/')[-1]
        ereg_hist2 = ereg_hist_filepath2.split('/')[-1]
        ereg_network2 = ereg_network_filepath2.split('/')[-1]
        partcoef_hist2 = partcoef_hist_filepath2.split('/')[-1]
        partcoef_network2 = partcoef_network_filepath2.split('/')[-1]
        moduleimage2 = module_filepath2.split('/')[-1]
        springimage2 = spring_filepath2.split('/')[-1]
        
    import time
    t = time.ctime()
    
    metrics_list = []
    for num in range(len(fullnames)):
        cur_metrics = []
        for count,region in enumerate(fullnames[num]):
            cur_metrics.append([region,
                                degs[num][count],\
                                ccs[num][count],\
                                bcs[num][count],\
                                partitions[num][count],\
                                e_regs[num][count],
                                part_coefs[num][count]]) #,e_locs[num][count]])
        metrics_list.append(cur_metrics)
    session.metrics_list = metrics_list

    d = dict(network1=network1,raw_density1=raw_densities[0],chosen_density1=chosen_densities[0],
                cmatimage1=cmatimage1,degreedistimage1=degreedistimage1,
                netimage1=netimage1,barimage1=barimage1,netimage2=netimage2,barimage2=barimage2,
                netimage3=netimage3,barimage3=barimage3,
                ereg_hist1=ereg_hist1,ereg_network1=ereg_network1,partcoef_hist1=partcoef_hist1,partcoef_network1=partcoef_network1,
                moduleimage1=moduleimage1,springimage1=springimage1,
                cpl1=cpl[0],mcc1=mcc[0],num_comp1=num_comp[0],
                sigma1=sigma[0],lambda_cc1=lambda_cc[0],gamma_cpl1=gamma_cpl[0],q1=q[0],
                network2=network2,raw_density2=raw_densities[1],chosen_density2=chosen_densities[1],
                cmatimage2=cmatimage2,degreedistimage2=degreedistimage2,
                netimage4=netimage4,barimage4=barimage4,netimage5=netimage5,barimage5=barimage5,
                netimage6=netimage6,barimage6=barimage6,
                ereg_hist2=ereg_hist2,ereg_network2=ereg_network2,partcoef_hist2=partcoef_hist2,partcoef_network2=partcoef_network2,
                moduleimage2=moduleimage2,springimage2=springimage2,
                cpl2=cpl[1],mcc2=mcc[1],num_comp2=num_comp[1],
                sigma2=sigma[1],lambda_cc2=lambda_cc[1],gamma_cpl2=gamma_cpl[1],q2=q[1], time=t,
                edge_attributes1=edge_attributes_list[0],edge_attributes2=edge_attributes_list[1])
    
    #report_filename = pdfreport(0, weight, network, chosen_density, raw_density, cpl,\
    #                        mcc, q, num_comp, sigma, lambda_cc, gamma_cpl,\
    #                        cmatimage, degreedistimage, barimage1, netimage1,\
    #                        barimage2, netimage2, barimage3, netimage3,\
    #                        moduleimage, springimage, metrics)
    #
    #d["report_filename"] = report_filename
    
    d["metrics_filename_1"] = metrics_text(metrics_list[0], network1, weights[0], chosen_densities[0])
    d["metrics_filename_2"] = metrics_text(metrics_list[1], network2, weights[1], chosen_densities[1])
    d["e_glob_1"] = e_glob[0]
    d["e_glob_2"] = e_glob[1] 
    
    row = db(db.analysis_results.network_name==network1.id).select().first()
    db(db.analysis_results.network_name==network1.id).update(num_times_analyzed=row.num_times_analyzed+1)
    now = datetime.datetime.now()
    db(db.analysis_results.network_name==network1.id).update(last_time_analyzed=now.strftime("%Y-%m-%d %H:%M"))
    row = db(db.analysis_results.network_name==network2.id).select().first()
    db(db.analysis_results.network_name==network2.id).update(num_times_analyzed=row.num_times_analyzed+1)
    now = datetime.datetime.now()
    db(db.analysis_results.network_name==network2.id).update(last_time_analyzed=now.strftime("%Y-%m-%d %H:%M"))
    
    return response.render(d)
  
def metrics():
    return dict(metrics=session.metrics, network=session.network)

def metrics_compare():
    return dict(metrics_list=session.metrics_list, network1=session.network1, network2=session.network2)
    
def metrics_text(metrics, network, weight, chosen_density):
    metrics_filename = os.path.join(network.network_name \
                 + '_' + weight + '_' + str(chosen_density) \
                 + '_' + 'metrics' + '.txt')
    metrics_filepath = os.path.join(request.folder,'static',metrics_filename)
    f = open(metrics_filepath, 'w')
    #f.write('#region_name\t#degree\t#clustering_coefficient\t#betweenness_centrality\t#module\t#regional_efficiency\t#local_efficiency\n')
    f.write('#region_name\t#degree\t#clustering_coefficient\t#betweenness_centrality\t#module\t#regional_efficiency\t#participation_coefficient\n')
    for region in metrics:
        s = '\t'.join([str(m) for m in region])
        s = s + '\n'
        f.write(s)
    f.close()
    return metrics_filename

# work in progres to create sortable metrics table using powertables
#def metrics_test():    
#    from gluon.sql import SQLTable
#    metrics = session.metrics
#    a=SQLTable('metrics') # fake table a
#    db.define_table('b',SQLField('name'),SQLField('a',a)) #reference it
#    db.define_table('metrics',SQLField('region_name'),SQLField('degree',db.b)) #define 
#    
#    #d = db.define_table(Field('region_name'),Field('degree'),Field('clustering_coefficient'),Field('betweenness_centrality'),Field('module'),Field('regional_efficiency'))
#    
#    #js_str = "'aa:data':["
#    #js_str = js_str + ','.join(str(x) for x in metrics)
#    #js_str = js_str + "],"
#    
#    #js_str_cols = "'aoColumns': [{ 'sTitle': 'Region Name' },{ 'sTitle': 'Degree' },{ 'sTitle': 'Clustering Coefficient' },{ 'sTitle': 'Betweenness Centrality' },{ 'sTitle': 'Module' },{ 'sTitle': 'Regional Efficiency' }]"
#    
#    #temp = {'region_name':'test'}
#    
#    #class Virtual(object):
#    #    @virtualsettings(label=T('region_name'))
#    #    def details_page(self):
#    #        return temp
#    
#    # could load in real table, hide all columns, and only have virtual; hack ^ n
#    
#    metricsTable = plugins.powerTable
#    metricsTable.datasource = a
#    #metricsTable.virtualfields = Virtual()
#    #metricsTable.dtfeatures['aaData'] = js_str
#    #metricsTable.dtfeatures['aoColumns'] = js_str_cols
#    #metricsTable.columns = ['virtual.region_name']
#
#    #print js_str
#    #print js_str_cols
#    return dict(table=metricsTable.create())

#def download():
#    import os
#    db = get_database(request)
#    return response.download(request,db)

def download():
    file = request.args[0]
    type = file.split('.')[1]
    if type == 'connectivity_matrix_file':
        network = session.network
        row = db(db.analysis_results.network_name==network.id).select().first()
        db(db.analysis_results.network_name==network.id).update(num_times_downloaded=row.num_times_downloaded+1)
        now = datetime.datetime.now()
        db(db.analysis_results.network_name==network.id).update(last_time_downloaded=now.strftime("%Y-%m-%d %H:%M"))
    return response.download(request, db)
    
def fastdownload():
    file = request.vars['file']
    return response.stream(open(file, 'rb'))

@auth.requires_login()
def upload():
    db.upload_data.email.default = auth.user.email
    #db.upload_data.email.writable = False
    form = SQLFORM(db.upload_data)
                   #keepopts=['study_name','network_name','email','atlas',
                   #                          'region_names_full_file',
                   #                          'region_names_abbrev_file',
                   #                          'region_xyz_centers_file',
                   #                          'connectivity_matrix_file',
                   #                          'imaging_modality','share','scanner_device','scan_parameters',
                   #                          'age_range_min','age_range_max','gender','subject_pool',
                   #                          'group_size','preprocessing_notes','funding_source'])
    form[0][0].insert(-1,"Succinct identifier for the location/purpose of the study, e.g. UCLA_ICBM")
    form[0][1].insert(-1,"Succinct name for inidividual matrix you're uploading, e.g. CONTROL_grpmean")
    form[0][2].insert(-1,"Your email address")
    form[0][3].insert(-1,"The Atlas/Method used to define regions, e.g. Freesurfer_68")
    form[0][4].insert(-1,"Text file listing full name of each region on separate lines")
    form[0][4].insert(-1, A('example', _href=URL('static', 'example_region_full_names.txt')))
    form[0][5].insert(-1,"Text file listing abbreviated name of each region on separate lines")
    form[0][5].insert(-1, A('example', _href=URL('static', 'example_region_abbrev_names.txt')))
    form[0][6].insert(-1,"Text file with (X,Y,Z) coordinate for each region")
    form[0][6].insert(-1, A('example', _href=URL('static', 'example_xyz_centers.txt')))
    form[0][7].insert(-1,"Text file with connectivity matrix")
    form[0][7].insert(-1, A('example', _href=URL('static', 'example_connectmat.txt')))
    form[0][8].insert(-1,"fMRI, DTI, DSI, MP-RAGE, EEG, MEG, ...")
    form[0][9].insert(-1,"Public data is viewable by anyone, Private is only viewable by you after logging in")
    form[0][10].insert(-1,"Siemens Trio 3T, GE Signa 3T, ...")
    form[0][11].insert(-1,"TR=2s, TE=30ms, Voxel Size=2x2x2, Diffusion directions=60, B value=1000, ...")
    form[0][15].insert(-1,"Alzheimer's, ASD, Schizophrenia, Normal, ...")
    form[0][16].insert(-1,"1 for an individual, otherwise the number of subjects averaged together")
    form[0][17].insert(-1,"Brief description of preprocessing steps/tools")
    form[0][18].insert(-1,"Funding source, publication reference, pubmed id, and relevant links")
    if form.accepts(request.vars, session):
        #f = request.vars.region_names_full_file.file.read()
        #print f
        #f1 = open(request.vars.region_names_full_file)
        #print f1.readlines()
        id = form.vars.id
        db.analysis_results.insert(network_name=id)
        response.flash = 'Your data has been uploaded'
    elif form.errors:
        response.flash = 'Form has errors'
    return dict(form=form)

@auth.requires_login()
def upload_batch():
    form = SQLFORM.factory(
        Field('study_name'),
        Field('network_names_list_file','upload', required=True, requires=IS_NOT_EMPTY(error_message='Select a file'), uploadfolder=os.path.join(request.folder,'uploads')),
        Field('email', required=True, default=auth.user.email),
        Field('atlas', required=True),
        Field('number_of_regions', 'integer', requires=IS_NOT_EMPTY()),
        Field('region_names_full_file_batch','upload', required=True, requires=IS_NOT_EMPTY(error_message='Select a file'), uploadfolder=os.path.join(request.folder,'uploads')),
        Field('region_names_abbrev_file_batch','upload', required=True, requires=IS_NOT_EMPTY(error_message='Select a file'), uploadfolder=os.path.join(request.folder,'uploads')),
        Field('region_xyz_centers_file_batch','upload', required=True, requires=IS_NOT_EMPTY(error_message='Select a file'), uploadfolder=os.path.join(request.folder,'uploads')),
        Field('connectivity_matrix_file_batch','upload', required=True, requires=IS_NOT_EMPTY(error_message='Select a file'), uploadfolder=os.path.join(request.folder,'uploads')),
        Field('imaging_modality', required=True, requires=IS_IN_SET(['DTI', 'DSI', 'HARDI', 'fMRI', 'Strucutral MRI', 'EEG', 'MEG', 'Other'])),
        Field('share', required=True, requires=IS_IN_SET(['Public', 'Private']), default='Public'),
        Field('scanner_device'),
        Field('scan_parameters'),
        Field('age_list_file','upload', uploadfolder=os.path.join(request.folder,'uploads')),
        Field('gender_list_file','upload', uploadfolder=os.path.join(request.folder,'uploads')),
        Field('subject_pool_list_file','upload', uploadfolder=os.path.join(request.folder,'uploads')),
        Field('group_size_list_file','upload', required=True, requires=IS_NOT_EMPTY(error_message='Select a file'), uploadfolder=os.path.join(request.folder,'uploads')),
        Field('preprocessing_notes','text'),
        Field('funding_source','text'))
    
    form[0][0].insert(-1,"Succinct identifier for the location/purpose of the study, e.g. UCLA_ICBM")
    form[0][1].insert(-1,"Text file listing succinct unique name for each network you're uploading, e.g. CONTROL_1, PATIENT_1, etc.")
    form[0][1].insert(-1, A('example', _href=URL('static', 'example_network_names_list.txt'), _title='Each network name should be listed on a separate line of the text file.'))
    form[0][2].insert(-1,"Your email address")
    form[0][3].insert(-1,"The Atlas/Parcellation Method used to define regions, e.g. Freesurfer_68")
    form[0][4].insert(-1,"The number of regions in the network")
    form[0][5].insert(-1,"Text file listing full name of each region on separate lines")
    region_names_text = "Each region name should be listed on a separate line. " \
    "If you're uploading six networks with 110 regions each, list all the region names " \
    "in order for the first network, then list all the region names again for the second " \
    "network, and so forth. In this example, the file will be 660 (6x110) lines " \
    "long. Make sure not to have any empty lines or extra characters between the " \
    "region names for the the different individual networks. In most cases, each " \
    "individual network will have identical region names. Click to see the example."
    form[0][5].insert(-1, A('example', _href=URL('static', 'example_region_full_names_batch.txt'), _title=region_names_text))
    form[0][6].insert(-1,"Text file listing abbreviated name of each region on separate lines")
    region_abbrevs_text = "Each region name abbreviation should be listed on a separate line. " \
    "If you're uploading six networks with 110 regions each, list all the region name abbreviations " \
    "in order for the first network, then list all the region name abbreviations again for the second " \
    "network, and so forth. In this example, the file will be 660 (6x110) lines " \
    "long. Make sure not to have any empty lines or extra characters between the " \
    "region name abbreviations for the the different individual networks. In most cases, each " \
    "individual network will have identical region name abbrevations. Click to see the example."
    form[0][6].insert(-1, A('example', _href=URL('static', 'example_region_abbrev_names_batch.txt'), _title=region_abbrevs_text))
    form[0][7].insert(-1,"Text file with (X,Y,Z) coordinate for each region")
    region_centers_text = "Each region xyz center should be space/tab delimited and listed on a separate line. " \
    "If you're uploading six networks with 110 regions each, list all the region xyz centers " \
    "in order for the first network, then list all the region xyz centers again for the second " \
    "network, and so forth. In this example, the file will be 660 (6x110) lines " \
    "long. Make sure not to have any empty lines or extra characters between the " \
    "region xyz for the the different individual networks. In some cases, each " \
    "individual network will have identical region xyz centers. In other cases, " \
    "different networks will have slightly different xyz centers based on individual " \
    "variations in structural or functional parcellation. PLEASE USE MNI mm COORDINATES IF POSSIBLE. Click to see the example."
    form[0][7].insert(-1, A('example', _href=URL('static', 'example_xyz_centers_batch.txt'), _title=region_centers_text))
    form[0][8].insert(-1,"Text file with connectivity matrix")
    connectmat_text = "Individual connectivity matrices should appear as space/tab " \
    "delimited text. Successive connectivity matrices should be stacked vertically " \
    "in the text file. For example, if you're uploading six connectivity matrices " \
    "with 110 regions each, the first connectivity matrix should appear on lines 1-110, " \
    "the second connectivity matrix on lines 111-220, and so forth. " \
    "In this example, the file will have 660 lines (ie rows) with 110 entries per line (ie columns). " \
    "Make sure not to have any empty lines or extra characters between the " \
    "connectivity matrices for the the different individual networks. In most cases, each " \
    "individual connectivity matrix will be different. Click to see the example."
    form[0][8].insert(-1, A('example', _href=URL('static', 'example_connectmat_batch.txt'), _title=connectmat_text))
    form[0][9].insert(-1,"fMRI, DTI, DSI, MP-RAGE, EEG, MEG, ...")
    form[0][10].insert(-1,"Public data is viewable by anyone, Private is only viewable by you after logging in")
    form[0][11].insert(-1,"Siemens Trio 3T, GE Signa 3T, ...")
    form[0][12].insert(-1,"TR=2s, TE=30ms, Voxel Size=2x2x2, Diffusion directions=60, B value=1000, ...")
    form[0][13].insert(-1,"Text file with age for each subject on seperate line")
    form[0][13].insert(-1, A('example', _href=URL('static', 'example_age_list.txt'), _title='The subject age corresponding to each individual network should be listed on a separate line of the text file.'))
    form[0][14].insert(-1,"Text file with gender for each subject on seperate line (accepts ONLY Male, Female, Mixed, or Unknown)")
    form[0][14].insert(-1, A('example', _href=URL('static', 'example_gender_list.txt'), _title='The gender age corresponding to each individual network should be listed on a separate line of the text file.'))
    form[0][15].insert(-1,"Text file with subject pool (Alzheimer's, ASD, Schizophrenia, Normal) for each subject on seperate line")
    form[0][15].insert(-1, A('example', _href=URL('static', 'example_subject_pool_list.txt'), _title='The subject status corresponding to each individual network should be listed on a separate line of the text file.'))
    form[0][16].insert(-1,"Text file with number of subjects used to calculate each individual network")
    groupsize_text = "The number of subjects used to calculate each individual network " \
                     "should be listed on a separate line of the text file. For networks " \
                     "calculated from individual subjects, this value should be 1. " \
                     "For networks calculated as an average of a group, this value should " \
                     "be the size of the group."
    form[0][16].insert(-1, A('example', _href=URL('static', 'example_group_size_list.txt'), _title=groupsize_text))
    form[0][17].insert(-1,"Brief description of preprocessing steps/tools")
    form[0][18].insert(-1,"Funding source, publication reference, pubmed id, and relevant links")

    if form.accepts(request.vars, session):
        num_regions = int(form.vars.number_of_regions)
        f1_path = request.folder + 'uploads/' + form.vars.network_names_list_file_newfilename
        f1 = open(f1_path)
        f1_items = [x.rstrip('\r\n') for x in f1]
        num_networks = len(f1_items)

        f2_path = request.folder + 'uploads/' + form.vars.region_names_full_file_batch_newfilename
        f2 = open(f2_path)
        f2_items = [x.rstrip('\r\n') for x in f2]
        f3_path = request.folder + 'uploads/' + form.vars.region_names_abbrev_file_batch_newfilename
        f3 = open(f3_path)
        f3_items = [x.rstrip('\r\n') for x in f3]
        f4_path = request.folder + 'uploads/' + form.vars.region_xyz_centers_file_batch_newfilename
        f4 = open(f4_path)
        f4_items = [x.rstrip('\r\n') for x in f4]
        f5_path = request.folder + 'uploads/' + form.vars.connectivity_matrix_file_batch_newfilename
        f5 = open(f5_path)
        f5_items = [x.rstrip('\r\n') for x in f5]
        
        f6_path = request.folder + 'uploads/' + form.vars.age_list_file_newfilename
        f6 = open(f6_path)
        f6_items = [x.rstrip('\r\n') for x in f6]
        f7_path = request.folder + 'uploads/' + form.vars.gender_list_file_newfilename
        f7 = open(f7_path)
        f7_items = [x.rstrip('\r\n') for x in f7]
        f8_path = request.folder + 'uploads/' + form.vars.subject_pool_list_file_newfilename
        f8 = open(f8_path)
        f8_items = [x.rstrip('\r\n') for x in f8]
        f9_path = request.folder + 'uploads/' + form.vars.group_size_list_file_newfilename
        f9 = open(f9_path)
        f9_items = [x.rstrip('\r\n') for x in f9]
        
        start = 0
        for i in range(num_networks):
            # split batch files into individual subject files
            filepath1 = request.folder + 'uploads/' + '%s_%s_region_names.txt' %(form.vars.study_name, f1_items[i])
            a = open(filepath1, 'w')
            for item in f2_items[start:start+num_regions]:
                a.write("%s\n" %item)
            a.close()
            stream1 = open(filepath1,'rb')
            
            filepath2 = request.folder + 'uploads/' + '%s_%s_region_abbrevs.txt' %(form.vars.study_name, f1_items[i])
            a = open(filepath2, 'w')
            for item in f3_items[start:start+num_regions]:
                a.write("%s\n" %item)
            a.close()
            stream2 = open(filepath2,'rb')
            
            filepath3 = request.folder + 'uploads/' + '%s_%s_region_xyz_centers.txt' %(form.vars.study_name, f1_items[i])
            a = open(filepath3, 'w')
            for item in f4_items[start:start+num_regions]:
                a.write("%s\n" %item)
            a.close()
            stream3 = open(filepath3,'rb')
            
            filepath4 = request.folder + 'uploads/' + '%s_%s_connectmat.txt' %(form.vars.study_name, f1_items[i])
            a = open(filepath4, 'w')
            for item in f5_items[start:start+num_regions]:
                a.write("%s\n" %item)
            a.close()
            stream4 = open(filepath4,'rb')
            
            start = start + num_regions
            
            try:
                id = db.upload_data.insert(study_name=form.vars.study_name,
                           network_name=f1_items[i],
                           email=form.vars.email,
                           atlas=form.vars.atlas,
                           region_names_full_file=db.upload_data.region_names_full_file.store(stream1,filepath1),
                           region_names_abbrev_file=db.upload_data.region_names_abbrev_file.store(stream2,filepath2),
                           region_xyz_centers_file=db.upload_data.region_xyz_centers_file.store(stream3,filepath3),
                           connectivity_matrix_file=db.upload_data.connectivity_matrix_file.store(stream4,filepath4),
                           imaging_modality=form.vars.imaging_modality,
                           share=form.vars.share,
                           scanner_device=form.vars.scanner_device,
                           scan_parameters=form.vars.scan_parameters,
                           age_range_min=f6_items[i],
                           age_range_max=f6_items[i],
                           gender=f7_items[i],
                           subject_pool=f8_items[i],
                           group_size=f9_items[i],
                           preprocessing_notes=form.vars.preprocessing_notes,
                           funding_source=form.vars.funding_source)
                db.analysis_results.insert(network_name=id)
            except:
                form.errors.network_names_list_file = 'Network names must each be unique for chosen study name'
                response.flash = 'Form has errors'
                return dict(form=form)
        response.flash = 'Your data has been uploaded'
    elif form.errors:
        response.flash = 'Form has errors'
    return dict(form=form)

def browse_OLD():
    """Browse all public data, allow users to update or delete data they've shared"""
    # can show results as all that are public, plus those that are created by the logged in user's email
    # add delete and update options to all that have user's email
    if auth.user:
        if request.vars:
            orderby = 'db' + '.' + str(request.vars.orderby)
            data = db((db.upload_data.share=='Public') | (db.upload_data.email==auth.user.email)).select(orderby=eval(orderby))
        else:
            data = db((db.upload_data.share=='Public') | (db.upload_data.email==auth.user.email)).select()
    else:
        data = db((db.upload_data.share=='Public')).select()
    return dict(data=data)

def update():
    """Allow user to update entries that they own, or view public entries in
    more detail"""
    row = db(db.upload_data.id == request.args[0]).select().first()
    session.network = row
    analysis_results_row = db(db.analysis_results.network_name==row.id).select().first()
    num_times_analyzed = analysis_results_row.num_times_analyzed
    last_time_analyzed = analysis_results_row.last_time_analyzed
    num_times_downloaded = analysis_results_row.num_times_downloaded
    last_time_downloaded = analysis_results_row.last_time_downloaded
    if auth.user:
        if row.email == auth.user.email:
            form = SQLFORM(db.upload_data, request.args[0], deletable=True,
                           delete_label=T('Check to delete'), url=URL('download'),\
                           submit_button='Update') 
        elif row.share == 'Public':
            form = SQLFORM(db.upload_data, request.args[0], readonly=True, upload=URL('download'))
    else:
        form = SQLFORM(db.upload_data, request.args[0], readonly=True, upload=URL('download')) ## EDIT url=URL..., this was a test
    email = row.email
    email_safe = email_make_safe(email)
    if form.accepts(request.vars, session):
        session.flash = T('done!')
        redirect(URL('browse'))
    return dict(form=form,
                num_times_analyzed=num_times_analyzed,
                last_time_analyzed=last_time_analyzed,
                num_times_downloaded=num_times_downloaded,
                last_time_downloaded=last_time_downloaded,
                email=email_safe)

def get_image(path, docheight, docwidth):
    """
    Set an image to fit on a page of a PDF, preserving aspect ratio
    """
    img = utils.ImageReader(path)
    width, height = img.getSize()
    aspect = height / float(width)
    #if height > docheight:
    if height > 700:
        #height = docheight
        height = 700
        width = height / aspect
    #elif width > docwidth:
    elif width > 500:
        #width = docwidth
        width = 500
        height = (width * aspect)
    return Image(path, width=width, height=height)

def pdfreport(num, weight, network, chosen_density, raw_density, cpl, mcc, q,\
              num_comp, sigma, lambda_cc, gamma_cpl,\
              cmatimage, degreedistimage, \
              barimage1, netimage1, \
              barimage2, netimage2, \
              barimage3, netimage3, \
              moduleimage, springimage, metrics):
    """
    Generate a PDF report containing Network Information, Global Network Metrics,
    Images, and Regional Network Metrics
    """
    from reportlab.rl_config import defaultPageSize
    PAGE_HEIGHT = defaultPageSize[1]; PAGE_WIDTH = defaultPageSize[0]
    
    styles = getSampleStyleSheet()
    
    report_filename = os.path.join(network.network_name \
                     + '_' + weight[num] + '_' + str(chosen_density) \
                     + '_' + 'report' + '.pdf')
    report_path = os.path.join(request.folder, 'static', report_filename)
    
    doc = SimpleDocTemplate(report_path,# pagesize=letter)
                            leftMargin=.5*inch,
                            rightMargin=.5*inch,
                            topMargin=.5*inch,
                            bottomMargin=.5*inch)
    report = []
    style = ParagraphStyle(
        name='Normal',
        fontSize=14
    )
    text = "<b>UCLA Multimodal Connectivity Database report</b>"
    p = Paragraph(text,style)
    report.append(p)
    report.append(Spacer(1,0.2*inch))
    
    headingstyle = ParagraphStyle(
        name='Normal',
        fontsize=12
    )
    
    # Network information
    text = "<b>Network Information</b>"
    style = styles["Normal"]
    p = Paragraph(text,style)
    report.append(p)
    
    text_items = []
    
    text="Network Name: %s" %network.network_name
    text_items.append(text)
    
    text="Atlas: %s" %network.atlas
    text_items.append(text)
    
    text="Imaging Modality: %s" %network.imaging_modality
    text_items.append(text)
    
    text="Subject Pool: %s" %network.subject_pool
    text_items.append(text)
    
    text="Group Size: %s" %network.group_size
    text_items.append(text)
    
    text="Age Range: %s-%s" %(network.age_range_min, network.age_range_max)
    text_items.append(text)
    
    text="Preprocessing Notes: %s" %network.preprocessing_notes
    text_items.append(text)
    
    for item in text_items:
        p = Paragraph(item, style)
        report.append(p)
        
    report.append(Spacer(1,0.2*inch))
    
    # Network global metrics
    text = "<b>Network Metrics</b>"
    p = Paragraph(text,style)
    report.append(p)
    
    gmetrics = []
    
    text="Raw Density: %s" %raw_density
    gmetrics.append(text)
    
    text="Chosen Density: %s" %chosen_density
    gmetrics.append(text)
    
    text="Characteristic Path Length: %s" %cpl
    gmetrics.append(text)
    
    text="Clustering Coefficient: %s" %mcc
    gmetrics.append(text)
    
    text="Modularity (Q): %s" %q
    gmetrics.append(text)
    
    text="Number of Components: %s" %num_comp
    gmetrics.append(text)
    
    text = "<i>Small World Attributes</i>"
    gmetrics.append(text)
    
    text="Sigma: %s" %sigma
    gmetrics.append(text)
    
    text="Lambda: %s" %lambda_cc
    gmetrics.append(text)
    
    text="Gamma: %s" %gamma_cpl
    gmetrics.append(text)
    
    for gmetric in gmetrics:
        p = Paragraph(gmetric, style)
        report.append(p)
    
    report.append(Spacer(1,0.2*inch))
    
    # Images
    text = "<b>Images</b>"
    p = Paragraph(text,style)
    report.append(p)
    
    imagelabelstyle = ParagraphStyle(
        name='Normal',
        fontsize=12,
        alignment=1)
    
    images = []
    image_texts = ["Connectivity Matrix", "Node Degree Distribution", \
                   "Node Degree Measures", "Network Degree Plot", \
                   "Node Betweenness Centrality Measures", "Network Betweenness Centrality Plot", \
                   "Node Clustering Measures", "Network Clustering Plot", \
                   "Network Modularity Plot", "Spring Modularity Plot"]
    image_files = [cmatimage, degreedistimage, \
                   barimage1, netimage1, \
                   barimage2, netimage2, \
                   barimage3, netimage3, \
                   moduleimage, springimage]

    for i in range(len(image_texts)):
        #img = Image(image_files[i]) #, width=300, height=300)
        #report.append(img)

        img = os.path.join(request.folder,'static', image_files[i])
        img = get_image(img, doc.height, doc.width)
        
        #report.append(get_image(image_files[i], doc.height, doc.width))
        report.append(img)
        
        p = Paragraph(image_texts[i], imagelabelstyle)
        report.append(p)
        report.append(Spacer(1,0.2*inch))
    
    # Regional Measures (create table)
    text = "<b>Regional Measures</b>"
    p = Paragraph(text,headingstyle)
    report.append(p)
    report.append(Spacer(1,0.2*inch))
    metrics_w_header = list(metrics)
    metrics_w_header.insert(0, ['Region Name','Degree','Clustering Coefficient','Betweenness Centrality','Module'])
    table = Table(metrics_w_header, [3*inch, .5*inch, 1.5*inch, 1.5*inch, .5*inch])
    table.setStyle(TableStyle([
    ('FONT', (0, 0), (-1, -1), 'Helvetica'),
    ('FONT', (0, 0), (-1, 0), 'Helvetica-Bold'),
    ('FONTSIZE', (0, 0), (-1, -1), 8),
    ('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black),
    ('BOX', (0, 0), (-1, 0), 0.25, colors.green),
    ('ALIGN', (0, 0), (-1, 0), 'CENTER')]))
    report.append(table)
    
    doc.build(report)
    return report_filename

def my_scoreatpercentile(in_mat, threshold_pct):
    '''
    Given a matrix and an integer percentage value (between 0-100),
    threshold the matrix to keep only the threshold_pct largest values,
    return the cutoff value
    '''
    a = in_mat.ravel()
    a_sort = a[a.argsort()[::-1]]
    dims = in_mat.shape
    dim_x = dims[0] # only need first dimension because matrices are always square
    if len(dims) > 1:
        thresh_value_pre = (threshold_pct/100.)*(dim_x*dim_x) # index of cutoff value
    else:
        thresh_value_pre = (threshold_pct/100.) * dim_x
    thresh_value = (round(thresh_value_pre / 2) * 2)
    if len(dims) > 1:
        if thresh_value >= dim_x*dim_x:
            cutoff = min(a_sort)
        else:
            cutoff = a_sort[thresh_value]
    else:
        if thresh_value >= dim_x:
            cutoff = min(a_sort)
        else:
            cutoff = a_sort[thresh_value]
    return cutoff

def email_make_safe(email):
    email_split = email.split('@')
    email_safe = email_split[0] + ' [AT] '
    email_split_post = email_split[-1].split('.')
    email_safe_post = ' [DOT] '.join(piece for piece in email_split_post)
    email_safe = email_safe + email_safe_post
    return email_safe

def browse():
    '''
    Browse all data either 1) publicly shared or 2) shared by the user who is logged in
    Uses web2py plugin powertables
    '''
    if auth.user:
        class Virtual(object):
            @virtualsettings(label=T('View/Download:'))
            def details_page(self):
                return A("View/Download", _href=URL("update", args=[self.upload_data.id]))
            @virtualsettings(label=T('View/Download:'))
            def analyze(self):
                return A("Analyze", _href=URL("index", vars={'study_name_cur':self.upload_data.study_name, 'network_name_cur':self.upload_data.network_name}))
            @virtualsettings(label=T('Information:'))
            def virtualtooltip(self):
                return T('Study: %s | Network: %s | Modality: %s' % (self.upload_data.study_name, self.upload_data.network_name, self.upload_data.imaging_modality))
            @virtualsettings(label=T('Shared by'))
            def shared_by(self):
                email = self.upload_data.email
                email_safe = email_make_safe(email)
                return email_safe
        
        powerTable = plugins.powerTable
        powerTable.datasource = data = db((db.upload_data.share=='Public') | (db.upload_data.email==auth.user.email)).select()
        powerTable.virtualfields = Virtual()
        powerTable.headers = 'fieldname:capitalize'
        powerTable.dtfeatures['bJQueryUI'] = True
        powerTable.uitheme = 'cupertino'
        powerTable.extrajs = dict(tooltip={'value':'vitualtooltip'})
        powerTable.dtfeatures['sScrollX'] = '100%'
        powerTable.dtfeatures['sScrollY'] = '100%'
        powerTable.dtfeatures['sPaginationType'] = request.vars.get('pager','full_numbers')
        powerTable.showkeycolumn = False
        powerTable.keycolumn = 'upload_data.id'
        powerTable.dtfeatures['iDisplayLength'] = 25
        powerTable.columns = ['virtual.details_page','virtual.analyze','upload_data.study_name','upload_data.network_name',
                           'upload_data.atlas','upload_data.imaging_modality','upload_data.share',
                           'upload_data.scanner_device','upload_data.age_range_min',
                           'upload_data.age_range_max','upload_data.gender','upload_data.subject_pool',
                           'upload_data.group_size','virtual.shared_by']
    else:
        class Virtual(object):
            @virtualsettings(label=T('View/Download:'))
            def details_page(self):
                return A("View/Download", _href=URL("update", args=[self.upload_data.id]))
            @virtualsettings(label=T('View/Download:'))
            def analyze(self):
                return A("Analyze", _href=URL("index", vars={'study_name_cur':self.upload_data.study_name, 'network_name_cur':self.upload_data.network_name}))
            @virtualsettings(label=T('Information:'))
            def virtualtooltip(self):
                return T('Study: %s | Network: %s | Modality: %s' % (self.upload_data.study_name, self.upload_data.network_name, self.upload_data.imaging_modality))
            @virtualsettings(label=T('Shared by'))
            def shared_by(self):
                email = self.upload_data.email
                email_safe = email_make_safe(email)
                return email_safe
        
        powerTable = plugins.powerTable
        powerTable.datasource = db((db.upload_data.share=='Public')).select()
        powerTable.virtualfields = Virtual()
        powerTable.headers = 'fieldname:capitalize'
        powerTable.dtfeatures['bJQueryUI'] = True
        powerTable.uitheme = 'cupertino'
        powerTable.extrajs = dict(tooltip={'value':'vitualtooltip'})
        powerTable.dtfeatures['sScrollX'] = '100%'
        powerTable.dtfeatures['sScrollY'] = '100%'
        powerTable.dtfeatures['sPaginationType'] = request.vars.get('pager','full_numbers')
        powerTable.showkeycolumn = False
        powerTable.keycolumn = 'upload_data.id'
        powerTable.dtfeatures['iDisplayLength'] = 25
        powerTable.columns = ['virtual.details_page','virtual.analyze','upload_data.study_name','upload_data.network_name',
                           'upload_data.atlas','upload_data.imaging_modality','upload_data.share',
                           'upload_data.scanner_device','upload_data.age_range_min',
                           'upload_data.age_range_max','upload_data.gender','upload_data.subject_pool',
                           'upload_data.group_size','virtual.shared_by']

    return dict(table=powerTable.create())

def browse_studies():
    if auth.user:
        studies = db((db.upload_data.id > 0) & (db.upload_data.email==auth.user.email)).select()
    else:
        studies = db((db.upload_data.id > 0) & (db.upload_data.share=='Public')).select()

    study_names = []
    
    for row in studies:
        study_names.append(row.study_name)

    study_names = list(set(study_names))

    studies = {}
    for study_name in study_names:
        if auth.user:
            study = db((db.upload_data.study_name==study_name) & (db.upload_data.email==auth.user.email)).select()
        else:
            study = db((db.upload_data.study_name==study_name) & (db.upload_data.share=='Public')).select()
        
        imaging_modalities = []
        subject_pools = []
        for row in study:
            imaging_modalities.append(row.imaging_modality)
            subject_pools.append(row.subject_pool)
        imaging_modalities = list(set(imaging_modalities))
        subject_pools = list(set(subject_pools))
        
        studies[study_name] = []
        studies[study_name].append(len(study))
        studies[study_name].append(imaging_modalities)
        studies[study_name].append(subject_pools)

    return dict(study_names=study_names, studies=studies)

def get_study_metadata():
    if len(request.args) > 0:
        study_name = request.args[0]
        study_name_nospaces = re.sub(' ', '_', study_name)
        if auth.user:
            study = db((db.upload_data.study_name==study_name) & (db.upload_data.email==auth.user.email)).select()
        else:
            study = db((db.upload_data.study_name==study_name) & (db.upload_data.share=='Public')).select()
        return dict(study=study)
    else:
        return dict(message='Please enter a study name')
        
def get_study_data():
    if len(request.args) > 0:
        study_name = request.args[0]
        study_name_nospaces = re.sub(' ', '_', study_name)
        if auth.user:
            study = db((db.upload_data.study_name==study_name) & (db.upload_data.email==auth.user.email)).select()
        else:
            study = db((db.upload_data.study_name==study_name) & (db.upload_data.share=='Public')).select()
        new_dir = os.path.join(request.folder,'uploads',study_name_nospaces)
        os.system('mkdir %s' %(new_dir))
        region_names_abbrev_files = []
        region_names_full_files = []
        region_xyz_centers_files = []
        connectivity_matrix_files = []
        for row in study:
            network_name = row.network_name
            
            file1 = os.path.join(request.folder, 'uploads', row.region_names_abbrev_file)
            new_region_names_abbrev_file = '%s_%s.txt' %(network_name, 'region_names_abbrev_file')
            cmd1 = 'cp %s %s/%s' %(file1, new_dir, new_region_names_abbrev_file)
            os.system(cmd1)

            file2 = os.path.join(request.folder, 'uploads', row.region_names_full_file)
            new_region_names_full_file = '%s_%s.txt' %(network_name, 'region_names_full_file')
            cmd2 = 'cp %s %s/%s' %(file2, new_dir, new_region_names_full_file)
            os.system(cmd2)
            
            file3 = os.path.join(request.folder, 'uploads', row.region_xyz_centers_file)
            new_region_xyz_centers_file = '%s_%s.txt' %(network_name, 'region_xyz_centers_file')
            cmd3 = 'cp %s %s/%s' %(file3, new_dir, new_region_xyz_centers_file)
            os.system(cmd3)
            
            file4 = os.path.join(request.folder, 'uploads', row.connectivity_matrix_file)
            new_connectivity_matrix_file = '%s_%s.txt' %(network_name, 'connectivity_matrix_file')
            cmd4 = 'cp %s %s/%s' %(file4, new_dir, new_connectivity_matrix_file)
            os.system(cmd4)
        zip_dir = os.path.join(request.folder,'uploads')
        zip_file = '%s_data.zip' % (study_name_nospaces)
        os.system('cd %s; zip -r %s %s' %(zip_dir, zip_file, study_name_nospaces))
        return response.stream(os.path.join(request.folder, 'uploads', zip_file))
    else:
        return dict(message='Please enter a study name in the format: http://umcd.humanconnectomeproject.org/get_study_data/<study_name>')