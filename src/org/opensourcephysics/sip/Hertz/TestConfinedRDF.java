package org.opensourcephysics.sip.Hertz;

import java.util.Random;

/**
 * Check ConfinedRDF using centers placed independently and uniformly
 * in the same box. For these random centers, g(r) should be near 1
 * wherever enough pairs have been counted.
 */
public class TestConfinedRDF {

    public static void main(String[] args) {

        // Each configuration contains N centers. We generate a completely
        // new set of N positions numberOfConfigurations times.
        int N = 108;
        int numberOfConfigurations = 2000;

        // The walls are at x = 0 and x = Lx.
        // The periodic y and z directions each have length side.
        double side = 33.077609611254374;
        double Lx = 20.077609611254374;

        double radBinWidth = 0.002;

        // Using the same seed generates the same random positions
        // each time this test is run.
        Random random = new Random(12345);

        double[] x = new double[N]; // 108 x-coordinate slots
        double[] y = new double[N]; // 108 y-coordinate slots
        double[] z = new double[N]; // 108 z-coordinate slots

        // Create the RDF calculator. It keeps references to the x, y, z
        // arrays above, so update() will read whatever positions are
        // currently stored in those arrays. No pairs are counted here
        ConfinedRDF rdf = new ConfinedRDF(
            x, y, z, side, Lx, radBinWidth, "_uniform_test_002"
        );

        // One pass through this outer loop creates and counts
        // one complete configuration of N centers.
        for (int configuration = 0; configuration < numberOfConfigurations; configuration++) {

            // Give each center a new, independent random position.
            for (int i = 0; i < N; i++) {
                // Scaling puts the center inside the box.
                x[i] = Lx * random.nextDouble();     // Between the walls
                y[i] = side * random.nextDouble();   // Periodic direction
                z[i] = side * random.nextDouble();   // Periodic direction
            }

            // All N centers now have positions. Count their pair
            // distances and add the results to the RDF histogram.
            // This also increases ConfinedRDF's count by one.
            rdf.update();
        }

        // Normalize the accumulated pair counts and write r and g(r).
        rdf.writeRDF();

        System.out.println(
            "Uniform-center RDF test complete: N=" + N
            + ", Lx=" + Lx
            + ", numberOfConfigurations=" + numberOfConfigurations
        );
    }
}