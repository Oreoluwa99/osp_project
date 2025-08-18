package org.opensourcephysics.sip.ch03;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.frames.*;

/**
 * LissajousApp calculates and plots a Lissajous figure.
 *
 * @author Alan Denton 03/17/10
 * @version 1.0
 */
public class LissajousApp extends AbstractSimulation {
  PlotFrame frame = new PlotFrame("x", "y", "Lissajous Figure");
  Lissajous lissajous = new Lissajous();

  /** 
   * Constructs a LissajousApp
   */
  public LissajousApp() {
    frame.addDrawable(lissajous);
//  frame.setPreferredMinMax(-1.1, 1.1, -1.1, 1.1);
//  frame.setSquareAspect(true);
  }

  /**
   * Initializes the simulation.
   */
  public void initialize() {
    lissajous.clear(); // clear display
    double ax = control.getDouble("ax"); // x-amplitude
    double ay = control.getDouble("ay"); // y-amplitude
    double xmax = 1.1*ax;
    double ymax = 1.1*ay;
    frame.setPreferredMinMax(-xmax, xmax, -ymax, ymax);
    double dt = control.getDouble("dt");
    double omegax = control.getDouble("omegax"); // angular frequency of x input
    double omegay = control.getDouble("omegay"); // angular frequency of y input
    double phix = control.getDouble("phix/pi");  // phase angle of x input
    double phiy = control.getDouble("phiy/pi");  // phase angle of y input
    phix = phix*Math.PI;  // convert phase angle to radians
    phiy = phiy*Math.PI;
// set instance variables:
    lissajous.dt = dt;
    lissajous.ax = ax; 
    lissajous.ay = ay; 
    lissajous.omegax = omegax; 
    lissajous.omegay = omegay; 
    lissajous.phix = phix; 
    lissajous.phiy = phiy; 
  }

  /**
   * Does a time step.
   */
  public void doStep() {
    lissajous.step();  // advance the state by one time step
    frame.setMessage("t = "+decimalFormat.format(lissajous.t));
/*
    double x = lissajous.x;
    double y = lissajous.y;
    double t = lissajous.t;
    control.println("t = "+t+" x = "+x+" y = "+y);
*/
  }

  /**
   * Resets the simulation.
   */
  public void reset() {
    control.setValue("ax",1);
    control.setValue("ay",1);
    control.setValue("omegax",2);
    control.setValue("omegay",3);
    control.setValue("phix/pi",0.);
    control.setValue("phiy/pi",0.);
//  control.setValue("phix/pi",1./6.);
//  control.setValue("phiy/pi",1./4.);
    control.setValue("dt", 0.05);
    enableStepsPerDisplay(true);
  }

  /**
   * Starts the Java application.
   * @param args  command line parameters
   */
  public static void main(String[] args) {
    SimulationControl.createApp(new LissajousApp());
  }
}
