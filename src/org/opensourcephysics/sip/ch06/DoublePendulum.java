package org.opensourcephysics.sip.ch06;
import java.awt.*;
import org.opensourcephysics.display.*;
import org.opensourcephysics.numerics.*;

/**
 * DoublePendulum models and animates a double pendulum.
 *
 * @author based on DampedDrivenPendulum by Wolfgang Christian, Jan Tobochnik, Harvey Gould
 * @version 1.1  revised 03/06/10 by Alan Denton
 */
public class DoublePendulum implements Drawable, ODE {
  double state[] = new double[5]; // [theta1, omega1, theta2, omega2, time]
  double g; // natural frequency squared
  int nstep;                      // dt = 2*PI/nstep
  Color color = Color.RED;
  int pixRadius = 6;
  RK4 odeMethod = new RK4(this);

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
   * The state array variables are: [theta, angular velocity, time].
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
    double theta1 = state[0];
    double omega1 = state[1];
    double theta2 = state[2];
    double omega2 = state[3];
    double sin1 = Math.sin(theta1);
    double sin2 = Math.sin(theta2);
    double cos1 = Math.cos(theta1);
    double cos2 = Math.cos(theta2);
    double sin12 = Math.sin(theta1-theta2);
    double cos12 = Math.cos(theta1-theta2);
    double cos12sq = Math.pow(cos12,2);

    rate[0] = state[1]; // state[1] = omega1
    rate[1] = (g*(sin2*cos12-2.*sin1)-(omega1*omega1*cos12+omega2*omega2)*sin12)/(2.-cos12sq);
    rate[2] = state[3]; // state[3] = omega2
    rate[3] = (2.*g*(sin1*cos12-sin2)+(omega2*omega2*cos12+2.*omega1*omega1)*sin12)/(2.-cos12sq);
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
