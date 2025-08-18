package org.opensourcephysics.sip.ch03;
import java.awt.*;
import org.opensourcephysics.display.*;
import org.opensourcephysics.numerics.*;

/**
 * Charge solves and displays the time evolution of a charged particle
 * in electric and magnetic fields (GTC, Problem 3.17)
 *
 * @author Alan Denton 02/10/10
 */
public class Charge implements Drawable, ODE {
  Circle circle = new Circle();
  Trail trail = new Trail();
  double[] state = new double[7]; // {x,vx,y,vy,b,e,t}
  EulerRichardson odeSolver = new EulerRichardson(this);

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
    trail.addPoint(state[0], state[2]); // x,y
  }

  /**
   * Sets the state.
   *
   * @param x
   * @param vx
   * @param y
   * @param vy
   * @param b  // magnetic field in z-direction
   * @param e  // electric field in y-direction
   */
  public void setState(double x, double vx, double y, double vy, double b, double e) {
    state[0] = x;
    state[1] = vx;
    state[2] = y;
    state[3] = vy;
    state[4] = b;
    state[5] = e;
    state[6] = 0;  // initial time
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
    double b = state[4]; 
    double e = state[5]; 
    rate[0] = vx;      // rate of change of x
    rate[1] = b*vy;    // rate of change of vx 
    rate[2] = vy;      // rate of change of y
    rate[3] = e-b*vx;  // rate of change of vy
    rate[4] = 0;       // constant B field
    rate[5] = 0;       // constant E field
    rate[6] = 1;       // dt/dt = 1
  }

  /**
   * Draws the charge and its path.  Required for Drawable interface.
   *
   * @param panel the drawing panel
   * @param g the graphics context
   */
  public void draw(DrawingPanel panel, Graphics g) {
    circle.setXY(state[0], state[2]);
    circle.draw(panel, g);
    trail.draw(panel, g);
  }
}
