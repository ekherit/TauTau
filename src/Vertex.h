/*
 * =====================================================================================
 *
 *       Filename:  Vertex.h
 *
 *    Description:  Adjust vertex of the tracks using db
 *
 *        Version:  1.0
 *        Created:  28.01.2019 16:50:02
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *   Organization:  Budker Insitute of Nuclear Physics
 *
 * =====================================================================================
 */

#ifndef IBN_BES_VERTEX_H
#define IBN_BES_VERTEX_H

#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/Bootstrap.h"

#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/Helix.h"

struct Vertex_t
{
  double rho,z,phi;
  Vertex_t()
  {
    rho = -9999;
    z   = -9999;
    phi = -9999;
  };

  Vertex_t (RecMdcTrack * track)
  {
    z  = track->z();
    phi = track->phi();
    rho = sqrt( track->x()*track->x()  +  track->y()*track->y() );
  };

  Vertex_t (RecMdcKalTrack * track)
  {
    z  = track->z();
    phi = track->phi();
    rho = sqrt( track->x()*track->x()  +  track->y()*track->y() );
  };

  void use_db(RecMdcTrack * mdcTrk)
  {
    //  Reconstruct the vertex 
    Hep3Vector xorigin(0,0,0);
    IVertexDbSvc*  vtxsvc;
    Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
    if(vtxsvc->isVertexValid())
    {
      double* dbv = vtxsvc->PrimaryVertex(); 
      double*  vv = vtxsvc->SigmaPrimaryVertex();  
      xorigin.setX(dbv[0]);
      xorigin.setY(dbv[1]);
      xorigin.setZ(dbv[2]);
    }
    // Vertex game. copy from rhophi analysis 
    double phi0=mdcTrk->helix(1);
    double xv=xorigin.x();
    double yv=xorigin.y();
    //double Rxy=(mdc.x[i]-xv)*cos(phi0)+(mdc.y[i]-yv)*sin(phi0);
    //mdc.r[i]=Rxy;
    HepVector a = mdcTrk->helix();
    HepSymMatrix Ea = mdcTrk->err();
    HepPoint3D point0(0.,0.,0.);   // the initial point for MDC recosntruction
    HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]); 
    VFHelix helixip(point0,a,Ea); 
    helixip.pivot(IP);
    HepVector vecipa = helixip.a();
    double  Rvxy0=fabs(vecipa[0]);  //the nearest distance to IP in xy plane
    double  Rvz0=vecipa[3];         //the nearest distance to IP in z direction
    double  Rvphi0=vecipa[1];
    rho=Rvxy0;
    z=Rvz0;
    phi=Rvphi0;
  }
};
#endif
