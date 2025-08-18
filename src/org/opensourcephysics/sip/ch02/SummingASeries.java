package org.opensourcephysics.sip.ch02;

public class SummingASeries{
    public static void main (String [] args){
    
	int N = 1000;
	double sum = 0;
	
	for (int i = 1; i < N; i++){
	    sum += (1.0/(i*i));
	}

	System.out.println(sum);
    }
}
