package org.opensourcephysics.sip.ch03;
import java.awt.*;
import org.opensourcephysics.display.*;
import org.opensourcephysics.numerics.*;

/**
 * CurveBall models the dynamics of a spinning projectile and forms a template for other simulations.
 *
 * @author Alan Denton 09/18/08
 */
public class CurveBall implements Drawable, ODE {
  static final double g = 9.8;
  double[] state = new double[12]; // {x,vx,wx,y,vy,wy,z,vz,wz,t,cd,cm}
  int pixRadius = 6;              // pixel radius for drawing of projectile
  EulerRichardson odeSolver = new EulerRichardson(this);

  public void setStepSize(double dt) {
    odeSolver.setStepSize(dt);
  }

  /**
   * Steps (advances) the time.
   *
   * @param dt the time step.
   */
  public void step() {
    odeSolver.step(); // do one time step using selected algorithm
  }

  /**
   * Sets the state.
   *
   * @param x
   * @param vx
   * @param wx
   * @param y
   * @param vy
   * @param wy
   * @param z
   * @param vz
   * @param wz
   */
  public void setState(double x, double vx, double wx, double y, double vy, double wy, double z, double vz, double wz, double cd, double cm) {
    state[0] = x;
    state[1] = vx;
    state[2] = wx;
    state[3] = y;
    state[4] = vy;
    state[5] = wy;
    state[6] = z;
    state[7] = vz;
    state[8] = wz;
    state[9] = 0;  // initial time
    state[10] = cd;
    state[11] = cm;
  }

  /**
   * Gets the state.  Required for ODE interface.
   * @return double[] the state
   */
  public double[] getState() {
    return state;
  }

  /**
   * Gets the rate.  Required for ODE interface
   * @param state double[] the state
   * @param rate double[]  the computed rate
   */
  public void getRate(double[] state, double[] rate) {
    double vx = state[1]; 
    double wx = state[2]; 
    double vy = state[4];
    double wy = state[5]; 
    double vz = state[7]; 
    double wz = state[8]; 
    double cd = state[10]; // drag coefficient
    double cm = state[11]; // Magnus coefficient
    double v = Math.sqrt(vx*vx+vy*vy+vz*vz); // speed
    rate[0] = vx; // rate of change of x
    rate[3] = vy; // rate of change of y
    rate[6] = vz; // rate of change of z
    rate[1] = -cd*v*vx+cm*(wy*vz-wz*vy);     // rate of change of vx 
    rate[4] = -cd*v*vy+cm*(wz*vx-wx*vz);     // rate of change of vx 
    rate[7] = -cd*v*vz+cm*(wx*vy-wy*vx)-g;   // rate of change of vy 
    rate[2] = 0;        // dwx/dt = 0
    rate[5] = 0;        // dwy/dt = 0
    rate[8] = 0;        // dwz/dt = 0
    rate[9] = 1;        //  dt/dt = 1
  }

  /**
   * Draws the projectile. Required for Drawable interface.
   *
   * @param drawingPanel
   * @param g
   */
  public void draw(DrawingPanel drawingPanel, Graphics g) {
    int xpix = drawingPanel.xToPix(state[0]);
    int ypix = drawingPanel.yToPix(state[3]);
//  int zpix = drawingPanel.zToPix(state[6]);
    g.setColor(Color.red);
    g.fillOval(xpix-pixRadius, ypix-pixRadius, 2*pixRadius, 2*pixRadius);
    g.setColor(Color.green);
    int xmin = drawingPanel.xToPix(-100);
    int xmax = drawingPanel.xToPix(100);
    int y0 = drawingPanel.yToPix(0);
    g.drawLine(xmin, y0, xmax, y0); // draw a line to represent the ground
  }
}
