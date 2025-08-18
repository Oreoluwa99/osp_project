package org.opensourcephysics.sip.ch06;
import org.opensourcephysics.controls.*;
import org.opensourcephysics.frames.*;

/**
 * DoublePendulumHamiltonApp solves Hamilton's equations for a double pendulum, 
   plots a phase space diagram and Poincare map, and draws an animation.
 *
 * @author based on PoincareApp by Wolfgang Christian, Jan Tobochnik, Harvey Gould 
 * @version 1.1  revised 04/01/20 by Alan Denton
 */
public class DoublePendulumHamiltonApp extends AbstractSimulation {
  final static double PI = Math.PI; // defined for brevity
  PlotFrame angle1 = new PlotFrame("time", "q1", "Theta 1");
  PlotFrame momentum1 = new PlotFrame("time", "p1", "Angular Momentum 1");
  PlotFrame totalEnergy= new PlotFrame("time", "energy", "Total Energy");
  PlotFrame phaseSpace = new PlotFrame("q1", "p1", "Phase Space Plot");
  PlotFrame poincare = new PlotFrame("q1", "p1", "Poincare Plot");
  DoublePendulumHamilton pendulum = new DoublePendulumHamilton();
  DisplayFrame displayFrame = new DisplayFrame("Pendulum");

  /**
   * Constructs the DoublePendulumHamiltonApp and intializes the display.
   */
  public DoublePendulumHamiltonApp() {
    phaseSpace.setMarkerShape(0, 6); // second argument indicates a pixel
    poincare.setMarkerSize(0, 1);    // smaller size gives better resolution
    poincare.setMarkerColor(0, java.awt.Color.RED);
    phaseSpace.setMessage("t = "+0);
    displayFrame.addDrawable(pendulum);
    displayFrame.setPreferredMinMax(-2.2, 2.2, -2.2, 2.2);
  }

  /**
   * Initializes the animation and clears the plots.
   */
  public void initialize() {
    double q1 = control.getDouble("theta1"); // initial angle 1
    double q2 = control.getDouble("theta2"); // initial angle 2          
    double energy = control.getDouble("Energy"); // initial energy
    double g = control.getDouble("g"); // acceleration due to gravity
    pendulum.g = g;
    double sin12sq = Math.pow(Math.sin(q1-q2), 2);
    double p1 = 0; // assume initial p1=0
    double p2 = Math.sqrt((energy-g*(3-2*Math.cos(q1)-Math.cos(q2)))*(1.+sin12sq)); // assume initial p1=0
    //double omega2 = p2/(1.-0.5*Math.pow(Math.cos(q1-q2), 2));
    //double omega1 = -0.5*omega2*Math.cos(q1-q2);
    pendulum.initializeState(new double[] {q1, p1, q2, p2, 0});
    pendulum.nstep = control.getInt("nstep"); 
    double dt = 2.*PI/pendulum.nstep; // step size
    pendulum.setStepSize(dt);
    clear();
  }

  /**
   * Clears the plots.
   */
  public void clear() {
    phaseSpace.clearData();
    poincare.clearData();
    phaseSpace.render();
    poincare.render();
  }

  /**
   * Does a step by advancing the time by time step dt.
   *
   */
  public void doStep() {
    double state[] = pendulum.getState();
    int nstep = pendulum.nstep;
    pendulum.step(); // advances the state by one time step
// Shift angles to the range [-PI,PI]:
    if(state[0] > PI) {
      state[0] = state[0]-2*PI;
    } else if(state[0] < -PI) {
      state[0] = state[0]+2*PI;
    }
    if(state[2] > PI) {
      state[2] = state[2]-2*PI;
    } else if(state[2] < -PI) {
      state[2] = state[2]+2*PI;
    }

    double q1 = state[0];
    double p1 = state[1];
    double q2 = state[2];
    double p2 = state[3];
    double time = state[4];
    double cos1 = Math.cos(q1);
    double cos2 = Math.cos(q2);
    double cos12 = Math.cos(q1-q2);
    double sin1 = Math.sin(q1);
    double sin2 = Math.sin(q2);
    double sin12 = Math.sin(q1-q2);
    double sin12sq = sin12*sin12;

    //double omega1 = (p1-cos12*p2)/(1+sin12sq);
    //double omega2 = (2*p2-cos12*p1)/(1+sin12sq);
    //double p1 = 2*omega1+omega2*cos12; // angular momentum
    //double p2 = omega2+omega1*cos12;
    double g = pendulum.g;
    double energy = (0.5*p1*p1+p2*p2-p1*p2*cos12)/(1+sin12sq)+g*(3-2*cos1-cos2); // energy

    angle1.append(0, time, q1);
    totalEnergy.append(0, time, energy);
    momentum1.append(0, time, p1);
    phaseSpace.append(0, q1, p1);

    double dt = 2*PI/pendulum.nstep; // step size

// Add a data point to Poincare map: 
    if(Math.abs(q2) < 0.01 && p2 > 0) {
      poincare.append(0, q1, p1);
    }
    phaseSpace.setMessage("t = "+decimalFormat.format(time));
    poincare.setMessage("t = "+decimalFormat.format(time));
    if(phaseSpace.isShowing()) {
      phaseSpace.render();
    }
    if(poincare.isShowing()) {
      poincare.render();
    }
  }

  /**
   * Resets all parameters to their defaults.
   */
  public void reset() {
    control.setValue("theta1", 0.);
    control.setValue("theta2", 0.);
    control.setValue("Energy", 15.);
    control.setValue("g", 9.8);
    control.setValue("nstep", 1000);
    enableStepsPerDisplay(true);
  }

  /**
   * Starts the Java application.
   * @param args  command line parameters
   */
  public static void main(String[] args) {
    SimulationControl control = SimulationControl.createApp(new DoublePendulumHamiltonApp());
    control.addButton("clear", "Clear");
  }
}
