package org.opensourcephysics.sip.ch06;
import java.awt.*;
import org.opensourcephysics.display.*;
import org.opensourcephysics.numerics.*;

/**
 * DoublePendulumHamilton models a double pendulum based on Hamilton's equations and animates the motion.
 *
 * @author based on DampedDrivenPendulum by Wolfgang Christian, Jan Tobochnik, Harvey Gould
 * @version 1.1  revised 04/01/20 by Alan Denton
 */
public class DoublePendulumHamilton implements Drawable, ODE {
  double state[] = new double[5]; // [q1, p1, q2, p2, time]
  double g; // natural frequency squared
  int nstep; // dt = 2*PI/nstep
  Color color = Color.RED;
  int pixRadius = 6;
  RK4 odeMethod = new RK4(this);
  //EulerRichardson odeMethod = new EulerRichardson(this);

  /**
   * Sets the ODESolver's time step.
   *
   * @param dt double
   */
  public void setStepSize(double dt) {
    odeMethod.setStepSize(dt);
  }

  /**
   * Steps (advances) the time.
   *
   * @param dt the time step.
   */
  public void step() {
    odeMethod.step(); // execute one RK4 step
  }

  /**
   * Initializes the state by copying the given array into the state array.
   * The state array variables are: [angles, angular momenta, time].
   *
   * @param newState double[]
   */
  void initializeState(double[] newState) {
    System.arraycopy(newState, 0, state, 0, 5);
  }

  /**
   * Gets the state.
   *
   * @return double[]
   */
  public double[] getState() {
    return state;
  }

  /**
   * Gets the rate using the given state.
   * @param state double[]
   * @param rate double[]
   */
  public void getRate(double state[], double rate[]) {
    double q1 = state[0];
    double p1 = state[1];
    double q2 = state[2];
    double p2 = state[3];
    double sin1 = Math.sin(q1);
    double sin2 = Math.sin(q2);
    double cos1 = Math.cos(q1);
    double cos2 = Math.cos(q2);
    double sin12 = Math.sin(q1-q2);
    double sin12sq = sin12*sin12;
    double cos12 = Math.cos(q1-q2);
    double K = (0.5*state[1]*state[1]+state[3]*state[3]-state[1]*state[3]*cos12)/(1+sin12sq); // kinetic energy

    rate[0] = (state[1]-cos12*state[3])/(1+sin12sq); // q1-dot
    rate[1] = -state[1]*state[3]*sin12/(1+sin12sq)+2*K*sin12*cos12/(1+sin12sq)-2*g*sin1; // p1-dot
    rate[2] = (2*state[3]-cos12*state[1])/(1+sin12sq); // q2-dot
    rate[3] = state[1]*state[3]*sin12/(1+sin12sq)-2*K*sin12*cos12/(1+sin12sq)-g*sin2; // p2-dot
    rate[4] = 1; // rate of change of time
  }

  /**
   * Draws the pendulum. Required for Drawable interface.
   *
   * @param drawingPanel
   * @param g
   */
  public void draw(DrawingPanel drawingPanel, Graphics g) {
    int xpivot = drawingPanel.xToPix(0);
    int ypivot = drawingPanel.yToPix(0);
    int x1pix = drawingPanel.xToPix(Math.sin(state[0]));
    int y1pix = drawingPanel.yToPix(-Math.cos(state[0]));
    int x2pix = drawingPanel.xToPix(Math.sin(state[0])+Math.sin(state[2]));
    int y2pix = drawingPanel.yToPix(-Math.cos(state[0])-Math.cos(state[2]));
    g.setColor(Color.black);
    g.drawLine(xpivot, ypivot, x1pix, y1pix);  // top rod
    g.drawLine(x1pix, y1pix, x2pix, y2pix);    // bottom rod 
    g.setColor(color);
    g.fillOval(x1pix-pixRadius, y1pix-pixRadius, 2*pixRadius, 2*pixRadius); // bobs
    g.fillOval(x2pix-pixRadius, y2pix-pixRadius, 2*pixRadius, 2*pixRadius); 
  }
}
