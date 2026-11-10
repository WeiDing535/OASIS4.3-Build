"""Geometry behind the 3D structure page.

The window draws the model itself on a Tk canvas, so this module only turns a
model into points and projects them: a Cα trace from the result preview (or all
atoms read from the file), one rotation matrix and an orthographic projection.
Keeping it free of Tk lets the maths be tested without a display.
"""
import math

import numpy as np

# One colour per chain, cycled for models with more chains than entries.
PALETTE=('#2f6f9f','#b5651d','#4b8b3b','#8e4b9e','#b03a48','#2e8b8b',
         '#9a7b1f','#5b5f9e','#a35a2a','#3f7a4f','#7a4f9e','#8a3f5f')


def chain_traces(preview):
    """Return ``[(chain, points)]`` of consecutive Cα positions.

    ``preview`` is the ``preview`` block of :func:`phasebuild.oasis_workflow.
    collect_results`; chains with fewer than two points cannot be drawn.
    """
    traces=[]
    for block in preview or []:
        rows=block.get('points') or []
        points=[]
        for row in rows:
            try:
                points.append([float(row[0]),float(row[1]),float(row[2])])
            except (TypeError,ValueError,IndexError):
                continue
        if len(points)<2:
            continue
        traces.append((str(block.get('chain') or ''),np.asarray(points,float)))
    return traces


def read_atoms(path):
    """Return ``[(chain, coordinates, elements)]`` for every atom of the model."""
    import gemmi

    structure=gemmi.read_structure(str(path))
    result=[]
    if not len(structure):
        return result
    for chain in structure[0]:
        coordinates=[]
        elements=[]
        for residue in chain:
            for atom in residue:
                coordinates.append((atom.pos.x,atom.pos.y,atom.pos.z))
                elements.append(atom.element.name)
        if coordinates:
            result.append((chain.name,np.asarray(coordinates,float),np.asarray(elements)))
    return result


def rotation(yaw,pitch):
    """Rotation matrix for a horizontal drag (yaw) and a vertical one (pitch)."""
    cy,sy=math.cos(yaw),math.sin(yaw)
    cp,sp=math.cos(pitch),math.sin(pitch)
    around_y=np.array([[cy,0.0,sy],[0.0,1.0,0.0],[-sy,0.0,cy]])
    around_x=np.array([[1.0,0.0,0.0],[0.0,cp,-sp],[0.0,sp,cp]])
    return around_x@around_y


def fit(points,width,height,margin=26):
    """Centre and scale (pixels per Å) that put every point inside the canvas."""
    points=np.asarray(points,float)
    if not len(points):
        return np.zeros(3),1.0
    center=points.mean(axis=0)
    radius=float(np.max(np.linalg.norm(points-center,axis=1)))
    radius=max(radius,1e-6)
    room=max(min(width,height)/2.0-margin,1.0)
    return center,room/radius


def project(points,matrix,center,scale,width,height):
    """Orthographic projection: x right, y up, z towards the viewer."""
    local=(np.asarray(points,float)-center)@matrix.T
    x=width/2.0+local[:,0]*scale
    y=height/2.0-local[:,1]*scale
    return x,y,local[:,2]
