package org.opensourcephysics.sip.ch03;
import java.awt.*;
import org.opensourcephysics.display3d.simple3d.*;
import org.opensourcephysics.frames.Display3DFrame;
import org.opensourcephysics.numerics.*;

/**
 * Charge3D solves and displays in a 3D box the time evolution of a 
 * charged particle in electric and magnetic fields (GTC, Problem 3.17)
 *
 * @author Alan Denton 02/10/10
 */
public class Charge3D implements ODE {
  double[] state = new double[13]; // {x,vx,y,vy,z,vz,ex,ey,ez,bx,by,bz,t}
  EulerRichardson odeSolver = new EulerRichardson(this);
  // Graphical elements
  Display3DFrame frame = new Display3DFrame("Charge in a box");
  ElementCircle ball;
  ElementTrail trail;
  double ballRadius = 0.1;
  double min = -2.0, max = 2.0;

  /**
   * Constructs a Charge3D
   */
  public Charge3D() {
    frame.setPreferredMinMax(min, max, min, max, min, max);
    ball = new ElementCircle();
    ball.setSizeXYZ(2*ballRadius, 2*ballRadius, 2*ballRadius);
    frame.addElement(ball);
    trail = new ElementTrail();
    trail.setMaximumPoints(100);
    trail.getStyle().setLineColor(java.awt.Color.RED);
    frame.addElement(trail);
  }

  public void setStepSize(double dt) {
    odeSolver.setStepSize(dt);
    trail.clear();
  }

  /**
   * Steps (advances) the time.
   *
   * @param dt the time step.
   */
  public void step() {
    odeSolver.step(); // do one time step using selected algorithm
    ball.setXYZ(state[0], state[2], state[4]); // x,y,z
    trail.addPoint(state[0], state[2], state[4]); // x,y,z
//  frame.setMessage("t = "+decimalFormat.format(state[12]));
  }

  /**
   * Sets the state.
   *
   * @param x
   * @param vx
   * @param y
   * @param vy
   * @param z
   * @param vz
   * @param ex
   * @param ey
   * @param ez
   * @param bx
   * @param by
   * @param bz
   */
  public void setState(double x, double vx, double y, double vy, double z, double vz, double ex, double ey, double ez, double bx, double by, double bz) {
    state[0] = x;
    state[1] = vx;
    state[2] = y;
    state[3] = vy;
    state[4] = z;
    state[5] = vz;
    state[6] = ex;
    state[7] = ey;
    state[8] = ez;
    state[9] = bx;
    state[10] = by;
    state[11] = bz;
    state[12] = 0;  // initial time
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
    double vy = state[3]; 
    double vz = state[5]; 
    double ex = state[6];
    double ey = state[7];
    double ez = state[8];
    double bx = state[9];
    double by = state[10];
    double bz = state[11];
    rate[0] = vx;      // rate of change of x
    rate[1] = ex+vy*bz-vz*by;    // rate of change of vx 
    rate[2] = vy;      // rate of change of y
    rate[3] = ey+vz*bx-vx*bz;  // rate of change of vy
    rate[4] = vz;       // constant B field
    rate[5] = ez+vx*by-vy*bx;  // rate of change of vz
    rate[6] = 0;       // constant E field
    rate[7] = 0;       
    rate[8] = 0;      
    rate[9] = 0;       // constant B field
    rate[10] = 0;       
    rate[11] = 0;      
    rate[12] = 1;       // dt/dt = 1
  }
}
