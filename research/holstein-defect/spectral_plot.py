"""Lossless vector heatmaps: merge neighboring cells with identical RGBA.

No interpolation, spectral averaging, rebinning or numerical-data changes.
The merged rectangles are exactly the colors of the original colormap cells.
"""
from __future__ import annotations
import re
from pathlib import Path
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib import colormaps


def edges(centers):
    c=np.asarray(centers);mid=(c[:-1]+c[1:])/2
    return np.r_[2*c[0]-mid[0],mid,2*c[-1]-mid[-1]]


def vector_map(ax,x,y,values,cmap,norm):
    a=np.asarray(values);xe=edges(x);ye=edges(y)
    if a.shape!=(len(x),len(y)) or not np.isfinite(a).all():raise ValueError('Invalid map')
    cm=colormaps[cmap];rgba=cm(norm(a));polygons=[];scalars=[]
    recovered=np.empty_like(rgba)
    for i,row in enumerate(rgba):
        cuts=np.r_[0,np.flatnonzero(np.any(row[1:]!=row[:-1],axis=1))+1,len(y)]
        for lo,hi in zip(cuts[:-1],cuts[1:]):
            polygons.append([(xe[i],ye[lo]),(xe[i+1],ye[lo]),(xe[i+1],ye[hi]),(xe[i],ye[hi])])
            scalars.append(a[i,lo]);recovered[i,lo:hi]=row[lo]
    if not np.array_equal(recovered,rgba):raise ValueError('Color-run merge changed the map')
    mesh=PolyCollection(polygons,cmap=cm,norm=norm,linewidths=0,edgecolors='none',antialiased=False,rasterized=False)
    mesh.set_array(np.asarray(scalars));ax.add_collection(mesh)
    mesh.vector_audit=dict(cells=int(a.size),rectangles=len(polygons),identical_rgba=True)
    return mesh


def save_figure(fig,stem,meshes=()):
    stem=Path(stem)
    fig.savefig(str(stem)+'.svg',metadata={'Date':None})
    svg=Path(str(stem)+'.svg')
    svg.write_text(re.sub(r'(<g id="(?:PolyCollection|QuadMesh)_[^"]+")>',r'\1 shape-rendering="crispEdges">',svg.read_text()))
    for mesh in meshes:mesh.set_edgecolor('face');mesh.set_linewidth(.15)
    fig.savefig(str(stem)+'.pdf',metadata={'CreationDate':None,'ModDate':None})
    for mesh in meshes:mesh.set_edgecolor('none');mesh.set_linewidth(0)
    fig.savefig(str(stem)+'.png',dpi=180)
