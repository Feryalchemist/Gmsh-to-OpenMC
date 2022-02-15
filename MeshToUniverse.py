'''
ALGORITHM
1. create all surface plane from mesh
2. create all
'''

import random
import math
import os
import openmc
import time
'''
-------------------------------------------
'''
import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default='browser'
'''
-------------------------------------------
'''

def plane_equation(a,b,c):
    v1 = [b[i]-a[i] for i in range(0,3)]
    v2 = [c[i]-a[i] for i in range(0,3)]
    n = [v1[1]*v2[2]-v1[2]*v2[1],
         v1[2]*v2[0]-v1[0]*v2[2],
         v1[0]*v2[1]-v1[1]*v2[0]]
    e = [n[0],n[1],n[2],
         (sum([a[i]*n[i] for i in range(0,3)]))]
    co = [(a[i]+b[i]+c[i])/3 for i in range(0,3)]
    return e,n,co

def tetrahedron(a,b,c,d):
    plane=[[b,c,d],[a,d,c],[a,b,d],[a,c,b]]
    cent = [sum([i[0] for i in [a,b,c,d]])/4,
            sum([i[1] for i in [a,b,c,d]])/4,
            sum([i[2] for i in [a,b,c,d]])/4]
    bound = [[max([i[0] for i in [a,b,c,d]]),min([i[0] for i in [a,b,c,d]])],
             [max([i[1] for i in [a,b,c,d]]),min([i[1] for i in [a,b,c,d]])],
             [max([i[2] for i in [a,b,c,d]]),min([i[2] for i in [a,b,c,d]])]]
    for n in bound:
        n += [abs(n[0]-n[1])]
    norm = []; equa = []; cnto = []
    for p in plane:
        e,n,co = plane_equation(p[0],p[1],p[2])
        cost = sum([(cent[i]-co[i])*n[i] for i in range(0,3)])
        if cost >=0:
            e+=[1]
        elif cost<0:
            e+=[-1]
        norm+=[n]
        equa+=[e]
        cnto+=[co]
    return equa,cent,bound

def draw_figure(norm,equa,cnto,cent,a,b,c,d):
    fig =go.Figure(data=[go.Mesh3d(x=[i[0] for i in [a,b,c,d]],
                                   y=[i[1] for i in [a,b,c,d]],
                                   z=[i[2] for i in [a,b,c,d]],
                                   i=[1,0,0,0],
                                   j=[2,3,1,2],
                                   k=[3,2,3,1],
                                   opacity = 0.5),
                         go.Scatter3d(x=[a[0],b[0],c[0],d[0],a[0],c[0],d[0],b[0]],
                                      y=[a[1],b[1],c[1],d[1],a[1],c[1],d[1],b[1]],
                                      z=[a[2],b[2],c[2],d[2],a[2],c[2],d[2],b[2]],
                                      name='Wireframe',
                                      mode='lines',
                                      line=dict(width = 5,dash='dot')),
                         go.Scatter3d(x=[cent[0]],
                                      y=[cent[1]],
                                      z=[cent[2]],
                                      name='centroid',
                                      mode='markers',
                                      marker=dict(color='rgb(255, 0,0)'))])
    for i in range(0,4):
        fig.add_trace(go.Scatter3d(
            x=[cnto[i][0],cnto[i][0]+0.1*norm[i][0]*-equa[i][4]],
            y=[cnto[i][1],cnto[i][1]+0.1*norm[i][1]*-equa[i][4]],
            z=[cnto[i][2],cnto[i][2]+0.1*norm[i][2]*-equa[i][4]],
            mode='lines',name='Fn_'+str(i+1),
            line=dict(width = 5,color='rgb('+str(255-i*50)+', '+str(i*50)+','+str(i*25)+')')))
        fig.add_trace(go.Cone(
            x=[cnto[i][0] + 0.98*(0.1*norm[i][0]*-equa[i][4])],
            y=[cnto[i][1] + 0.98*(0.1*norm[i][1]*-equa[i][4])],
            z=[cnto[i][2] + 0.98*(0.1*norm[i][2]*-equa[i][4])],
            u=[0.5*(0.1*norm[i][0]*-equa[i][4])],
            v=[0.5*(0.1*norm[i][1]*-equa[i][4])],
            w=[0.5*(0.1*norm[i][2]*-equa[i][4])],
            colorscale=[[0, 'rgb('+str(255-i*50)+', '+str(i*50)+','+str(i*25)+')'], 
                        [1, 'rgb('+str(255-i*50)+', '+str(i*50)+','+str(i*25)+')']],
            showlegend=False,
            showscale=False
            ))
    fig.show()

def read_msh(filename,f=1):
    nodes = []; index_nodes = []; elements=[];
    with open(filename,'r') as file:
        nodes_cnt=False; elements_cnt=0
        for num,line in enumerate(file, 1):
            if '$Nodes\n' in line:
                nodes_cnt=True
            elif '$EndNodes\n' in line:
                nodes_cnt=False
            elif '$Nodes\n' not in line and nodes_cnt==True:
                if len(line.split())==3:
                    nodes.append([float(i)*f for i in line.split()])
                elif len(line.split())==1:
                    index_nodes+=[int(line)]
            if '$Elements\n' in line:
                elements_cnt = num+2
            if elements_cnt!=0:
                if elements_cnt==num:
                    if '$EndElements' not in line:
                        elements_cnt+=int(line.split()[-1])+1
                        elements += [[]]
                elif elements_cnt!=num and len(elements)!=0:
                    elements[-1] += [[int(i) for i in line.split()]]
    nodes = [[index_nodes[i]]+nodes[i] for i in range(0,len(nodes))]
    return nodes,elements

def surface_mesh(r,nodes,elements):
    surface_equation=[]
    for tri in elements:
        a,b,c = [nodes[[i[0] for i in nodes].index(tri[1])][1:],
                 nodes[[i[0] for i in nodes].index(tri[2])][1:],
                 nodes[[i[0] for i in nodes].index(tri[3])][1:]]
    
        e,n,co=plane_equation(a,b,c)
        th =  math.acos(co[0]/(co[0]**2+co[1]**2)**0.5)
        vr = [r*math.cos(th)-co[0],r*math.sin(th)-co[1],-co[2]]
        if sum([vr[i]*n[i] for i in range(0,3)]) >=1:
            e+=[-1]
        else:
            e+=[1]
        surface_equation+=[e]
    
    surface_equation = {'number':[t[0] for t in elements],'equation':surface_equation,
                        'nodes':[[t[1],t[2],t[3]] for t in elements]}
    return surface_equation

# This one will create collection of each tetrahedron surface from .msh
def fill_mesh(elements,nodes,surface_nodes={'number':[],'equation':[],'nodes':[[],[],[]]}):
    mesh_equation = []; bound_surface=[]; nodess=[]; cnt=0; centroid=[]; bbox=[]
    for q in elements:
        a,b,c,d = [nodes[[i[0] for i in nodes].index(q[1])],
                   nodes[[i[0] for i in nodes].index(q[2])],
                   nodes[[i[0] for i in nodes].index(q[3])],
                   nodes[[i[0] for i in nodes].index(q[4])]]
        equa,cent,bound= tetrahedron(a[1:],b[1:],c[1:],d[1:])
        mesh_equation += [equa]
        centroid += [cent]
        bbox += [bound]
        bs=[0]*4
        for num,j in enumerate([[b[0],c[0],d[0]],[a[0],d[0],c[0]],[a[0],b[0],d[0]],[a[0],c[0],b[0]]]):
            for k in surface_nodes:
                if set(j)==set(k):
                    bs[num] = 1; cnt+=1
                    #print(set(k),'---',set(j))
        bound_surface+=[bs]
        nodess+=[[[b[0],c[0],d[0]],[a[0],d[0],c[0]],[a[0],b[0],d[0]],[a[0],c[0],b[0]]]]
    mesh_equation = {'number':[q[0] for q in elements],'equation':mesh_equation,
                     'nodes':nodess,'bound_surface':bound_surface,
                     'bound_box':bbox,'centroid':centroid}
    return mesh_equation

# Mode 2d will create single cell from 2D surface collection from surface_mesh
# Mode 3D will create multiple cell from tetrahedron collection from filled_mesh
#         if there are any tetrahedron surface samilar to any 2D surface collection, then
#         it will be vacuum boundary

def create_surface(elements,mode='2d'):
    surface=[]
    for num, plane in enumerate(elements,1):
        surface += openmc.Plane(a=plane[0],b=plane[2],c=plane[3],d=plane[4])
        if plane[4]==1:
            b = +surface[num-1]
        else: 
            b = +surface[num-1]
        if num==1:
            bound = b
        else:
            bound = bound & b
    
    return surface,bound
            
def Plot_iMesh(fills_mesh):
    Plasma = openmc.Material(name = 'Plasma')
    Plasma.set_density('atom/b-cm',1e-5)
    Plasma.add_nuclide('H2',0.5)
    Plasma.add_nuclide('H3',0.5)
    material = openmc.Materials([Plasma])
    material.export_to_xml()

    NTetra = random.choice(list(range(0,len(fills_mesh['equation']))))
    S  = [];
    for i,eq in enumerate(fills_mesh['equation'][NTetra],0):
        S += [openmc.Plane(a=eq[0],b=eq[1],c=eq[2],d=eq[3],boundary_type='vacuum')]
        if i==0:
            if eq[4]==-1:
                region = -S[i];
            else:
                region = S[i];    
        else:
            if eq[4]==-1:
                region = region & -S[i]
            else:
                region = region & S[i]

    root_cell = openmc.Cell(name='Plasma',fill=Plasma,region=region)
    geometry = openmc.Geometry()
    geometry.root_universe= openmc.Universe(cells=[root_cell])
    geometry.export_to_xml()

    root_cell = openmc.Cell(name='Plasma',fill=Plasma,region=region)
    geometry = openmc.Geometry()
    geometry.root_universe= openmc.Universe(cells=[root_cell])
    geometry.export_to_xml()

    settings = openmc.Settings()
    settings.batches = 550
    settings.inactive = 10
    settings.particles = int(100000)
    settings.export_to_xml()

    res = 500

    plot1 = openmc.Plot(plot_id=1)
    plot1.filename = 'HTR10-core-xz'
    plot1.basis='xz'
    plot1.origin = fills_mesh['centroid'][NTetra]
    plot1.width = [fills_mesh['bound_box'][NTetra][0][2],fills_mesh['bound_box'][NTetra][2][2]]
    plot1.pixels = [res,int(res/fills_mesh['bound_box'][NTetra][0][2]*fills_mesh['bound_box'][NTetra][2][2])]
    plot1.color_by = 'material'

    plot2 = openmc.Plot(plot_id=2)
    plot2.filename = 'HTR10-core-xy'
    plot2.basis='xy'
    plot2.origin = fills_mesh['centroid'][NTetra]
    plot2.width = [fills_mesh['bound_box'][NTetra][0][2],fills_mesh['bound_box'][NTetra][1][2]]
    plot2.pixels = [res,int(res/fills_mesh['bound_box'][NTetra][0][2]*fills_mesh['bound_box'][NTetra][1][2])]
    plot2.color_by = 'material'

    openmc.plot_inline(plot1)
    openmc.plot_inline(plot2)
    
class Create_Cell_set():
    def __init__(self,elements_dict,nodes_dict,material):
        self.material = material
        self.mesh_dict= fill_mesh(elements_dict,nodes_dict)
    def Create(self):
        Surf = []; Cell = []
        for j,c in enumerate(self.mesh_dict['equation']):
            S = []
            for i,eq in enumerate(c,0):
                S += [openmc.Plane(a=eq[0],b=eq[1],c=eq[2],d=eq[3])]
                if i==0:
                    if eq[4]==-1:
                        region = -S[i];
                    else:
                        region = S[i];    
                else:
                    if eq[4]==-1:
                        region = region & -S[i]
                    else:
                        region = region & S[i]
            Surf += [S]
            Cell += [openmc.Cell(name='Plasma',fill=self.material,region=region)]
        self.mesh_dict['surface'] = Surf
        self.mesh_dict['cell']    = Cell

#Plot_iMesh(fills_mesh)
'''
Plasma = openmc.Material(name = 'Plasma')
Plasma.set_density('atom/b-cm',1e-5)
Plasma.add_nuclide('H2',0.5)
Plasma.add_nuclide('H3',0.5)
material = openmc.Materials([Plasma])
material.export_to_xml()

Surf  = []; Cell = []; iRegion=0
for j,c in enumerate(fills_mesh['equation']):
    S = []
    for i,eq in enumerate(c,0):
        # If it is a boundary surface then assign the region
        if fills_mesh['bound_surface'][j][i]==1:
            S += [openmc.Plane(a=eq[0],b=eq[1],c=eq[2],d=eq[3],boundary_type='vacuum')]
            if iRegion==0:
                iRegion = -S[i]; oRegion = +S[i]
            else:
                iRegion = iRegion & -S[i]; oRegion = oRegion | +S[i]
        else:
            S += [openmc.Plane(a=eq[0],b=eq[1],c=eq[2],d=eq[3])]
        # Checking the normal vector of each side of Tetrahedron
        if i==0:
            if eq[4]==-1:
                region = -S[i];
            else:
                region = S[i];    
        else:
            if eq[4]==-1:
                region = region & -S[i]
            else:
                region = region & S[i]
    Surf += S
    Cell += [openmc.Cell(name='Plasma',fill=Plasma,region=region)]
'''

#==============================================================================



#surface,bound = create_surface(elements)

#os.environ['OPENMC_CROSS_SECTIONS']='$HOME/Desktop/Nuclear_Library/jeff-3.3'
#plasma = openmc.Material()
#plasma.add_nuclide('H2',0.5,'ao')
#plasma.add_nuclide('H3',0.5,'ao')
#materials_file=openmc.Materials([plasma])

#Stellarator = openmc.Universe()
#Plasma = openmc.Cell(name='plasma',fill=plasma,region=bound)
#Stellarator.add_cell

    

    
    
    
        
        
        
    