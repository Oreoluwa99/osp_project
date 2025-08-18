package org.opensourcephysics.sip.ch04;
import org.opensourcephysics.numerics.ODE;
import java.awt.*;
import org.opensourcephysics.display.*;
import org.opensourcephysics.numerics.*;

/**
 * DDO implements the ODE interface for a driven, damped oscillator.
 *
 * @author Alan Denton
 * @version 1.0
 */
public class DDO implements Drawable, ODE {
  double x0, v0, omega0, omega0Squared, gamma, omega, a0, energy; 
  double[] state = new double[3]; // [displacement, velocity, time] state array
  Color color = Color.RED;
  int pixRadius = 6;

  /**
   * Construct a driven damped oscillator with the given values.
   */
  public DDO(double _x0, double _v0, double _omega0, double _gamma, double _omega, double _a0) {
    x0 = _x0;
    v0 = _v0;
    omega0 = _omega0;
    omega0Squared = omega0*omega0;
    gamma = _gamma;
    omega = _omega;
    a0 = _a0;
    state[0] = x0; // displacement
    state[1] = v0; // velocity
    state[2] = 0; // time
  }

  /**
   * Get state array. Implementation of ODE interface.
   *
   * @return state array
   */
  public double[] getState() {
    return state;
  }

  /**
   * Get the driving force at given time.
   * @param t
   * @return
   */
  public double getDrivingForce(double t) {
    return a0*Math.cos(omega*t);
  }

  /**
   * Get the rate array. Implementation of ODE interface.
   * This method may be invoked many times with different intermediate states
   * as an ODESolver is carrying out the solution.
   *
   * @param state the state array
   * @param rate the rate array
   */
  public void getRate(double[] state, double[] rate) {
    omega0Squared = omega0*omega0;
    rate[0] = state[1];
    double drivingForce = getDrivingForce(state[2]);  
    rate[1] = -omega0Squared*state[0]-gamma*state[1]+drivingForce; 
    rate[2] = 1;  
    double x = state[0];
    double v = state[1];
    energy = (v*v+omega0Squared*x*x)/2;
  }

  /**
   * Draws the oscillator. Required for Drawable interface.
   *
   * @param drawingPanel
   * @param g
   */
  public void draw(DrawingPanel drawingPanel, Graphics g) {
    int xanchor = drawingPanel.xToPix(0);
    int yanchor = drawingPanel.yToPix(0);
    int xpix = drawingPanel.xToPix(state[0]);
    int ypix = yanchor;
    g.setColor(Color.black);
    g.drawLine(xanchor, yanchor, xpix, ypix);                // the spring
    g.setColor(color);
    g.fillOval(xpix-pixRadius, ypix-pixRadius, 2*pixRadius, 2*pixRadius); // bob
  }
}
