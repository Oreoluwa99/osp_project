package org.opensourcephysics.sip.ch03;
import java.awt.*;
import org.opensourcephysics.display.*;

/**
 * Lissajous calculates and plots a Lissajous figure.
 *
 * @author Alan Denton 03/17/10
 * @version 1.0
 */
public class Lissajous implements Drawable {
  double x;
  double y;
  double omegax;
  double omegay;
  double phix;
  double phiy;
  double ax;
  double ay;
  double t = 0.;
  double dt;
  Trail trail = new Trail();

  /**
   * Steps (advances) the time.
   */
  public void step() {
    x = ax*Math.sin(omegax*t+phix);
    y = ay*Math.sin(omegay*t+phiy);
    t += dt;
    trail.addPoint(x, y); // x,y
  }

  /**
   * Clears the display.
   */
  public void clear() {
    trail.clear();
  }

  /**
   * Draws the Lissajous path.  Required for Drawable interface.
   *
   * @param panel the drawing panel
   * @param g the graphics context
   */
  public void draw(DrawingPanel panel, Graphics g) {
    trail.draw(panel, g);
  }
}
