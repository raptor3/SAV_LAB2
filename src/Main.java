import java.util.ArrayList;
import java.util.Collection;

import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.linear.LinearConstraintSet;
import org.apache.commons.math3.optim.linear.LinearObjectiveFunction;
import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.Relationship;
import org.apache.commons.math3.optim.linear.SimplexSolver;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;



public class Main {
	public static void main(String[] args) {
		final FuzzyNum __1 = new TrupFuzzyNum(1, 1, 1, 1);
		final FuzzyNum _1 = new TrupFuzzyNum(1, 1, 1, 3);
		final FuzzyNum _2 = new TrupFuzzyNum(1, 2, 2, 4);
		final FuzzyNum _3 = new TrupFuzzyNum(1, 3, 3, 5);
		final FuzzyNum _4 = new TrupFuzzyNum(2, 4, 4, 6);
		final FuzzyNum _5 = new TrupFuzzyNum(3, 5, 5, 7);
		final FuzzyNum _6 = new TrupFuzzyNum(4, 6, 6, 8);
		final FuzzyNum _7 = new TrupFuzzyNum(5, 7, 7, 9);
		final FuzzyNum _8 = new TrupFuzzyNum(6, 8, 8, 9);
		final FuzzyNum _9 = new TrupFuzzyNum(7, 9, 9, 9);		
		final FuzzyNum _1_1 = new InverseTrupFuzzyNum(1, 1, 1, 3);
		final FuzzyNum _1_2 = new InverseTrupFuzzyNum(1, 2, 2, 4);
		final FuzzyNum _1_3 = new InverseTrupFuzzyNum(1, 3, 3, 5);
		final FuzzyNum _1_4 = new InverseTrupFuzzyNum(2, 4, 4, 6);
		final FuzzyNum _1_5 = new InverseTrupFuzzyNum(3, 5, 5, 7);
		final FuzzyNum _1_6 = new InverseTrupFuzzyNum(4, 6, 6, 8);
		final FuzzyNum _1_7 = new InverseTrupFuzzyNum(5, 7, 7, 9);
		final FuzzyNum _1_8 = new InverseTrupFuzzyNum(6, 8, 8, 9);
		final FuzzyNum _1_9 = new InverseTrupFuzzyNum(7, 9, 9, 9);
		
        FuzzyNum[][] dc = new FuzzyNum[][] {
        		{ __1 , _1 , _5 },
        		{ _1_1, __1, _3 },
        		{ _1_5,_1_3, __1}
        };
        FuzzyNum[][] d1 = new FuzzyNum[][] {
        		{ 	__1, _1_7, _1_3, _1_2 },
        		{ 	_7, __1,	_3,	_3 	},
        		{ 	_3,	_1_3, __1, 	_1	},
        		{ 	_2,	_1_3, _1_1,	__1	}
        };
        FuzzyNum[][] d2 = new FuzzyNum[][] {
        		{ 	__1, _3, _5, _7 },
        		{ 	_1_3, __1,	_3,	_3 	},
        		{ 	_1_5,	_1_3, __1, 	_1	},
        		{ 	_1_7,	_1_3, _1_1,	__1	}
        };
        FuzzyNum[][] d3 = new FuzzyNum[][] {
        		{ 	__1, 	_2, _1_3,	 _3 },
        		{ 	_1_2, __1,	_1_3,	_1 	},
        		{ 	_3,		_3, __1, 	_2	},
        		{ 	_1_3,	_1_1, _1_2,	__1	}
        };
        
        Matrix Dc0 = new Matrix("Dc", dc, 0);
        Dc0.findW();
        Matrix D10 = new Matrix("D1", d1, 0);
        D10.findW();
        Matrix D20 = new Matrix("D2", d2, 0);
        D20.findW();
        Matrix D30 = new Matrix("D3", d3, 0);
        D30.findW();
        Matrix Dc05 = new Matrix("Dc", dc, 0.5);
        Dc05.findW();
        Matrix D105 = new Matrix("D1", d1, 0.5);
        D105.findW();
        Matrix D205 = new Matrix("D2", d2, 0.5);
        D205.findW();
        Matrix D305 = new Matrix("D3", d3, 0.5);
        D305.findW();
	}
}

interface FuzzyNum {
	
	public Interval getAlphaInterval(double alpha);
}

class TrupFuzzyNum implements FuzzyNum{
	double l;
	double m1;
	double m2;
	double u;
	
	public TrupFuzzyNum(double l, double m1, double m2, double u) {
		this.l = l;
		this.m1 = m1;
		this.m2 = m2;
		this.u = u;
	}
	
	public Interval getAlphaInterval(double alpha) {
		return new Interval(l+alpha*(m1-l), u-alpha*(u-m2));
	}
}

class InverseTrupFuzzyNum implements FuzzyNum{
	double l;
	double m1;
	double m2;
	double u;
	
	public InverseTrupFuzzyNum(double l, double m1, double m2, double u) {
		this.l = l;
		this.m1 = m1;
		this.m2 = m2;
		this.u = u;
	}
	
	public Interval getAlphaInterval(double alpha) {
		return new Interval(1/(u-alpha*(u-m2)),1/(l+alpha*(m1-l)));
	}
}

class Interval {
	public double L;
	public double U;
	
	public Interval(double l, double u) {
		L = l;
		U = u;
	}
	
	public String toString() {
		return "(" + String.format("%.3f", L)  + ";" + String.format("%.3f", U) + ")";
	}
}

class Matrix {
	public Interval[][] mass;
	private String name;
	private double a;
	
	public Matrix (String name, FuzzyNum[][] fuzzyNumb, double alpha) {
		this.name = name + " alpha " + alpha;
		mass = new Interval[fuzzyNumb.length][fuzzyNumb.length];
		for (int i = 0; i < mass.length; i++) {
			for (int j = 0; j < mass.length; j++) {
				mass[i][j] = fuzzyNumb[i][j].getAlphaInterval(alpha);
			}
		}
		System.out.println(toString());
	}
	
	public boolean MatrixIsCorrect() {
		boolean correct = true;
		for (int i = 0; i < mass.length; i++ ) {
			for (int j = 0; j < mass.length; j++) {
				double maxL = 0;
				double minU = Double.POSITIVE_INFINITY;
				for (int k =0; k < mass.length; k++) {
					if (minU > mass[i][k].U * mass[k][j].U) {
						minU = mass[i][k].U * mass[k][j].U;
					}
					if (maxL < mass[i][k].L * mass[k][j].L) {
						maxL = mass[i][k].L * mass[k][j].L;
					}
				}
				correct = correct && (maxL<=minU);
			}				
		}
		System.out.println("Matrix correct " + correct);
		if (!correct) {
        	FindCorrectMatrix();
        	System.out.println("correct one");
            System.out.println(toString());
        }		
		return correct;
	}
	
	public void FindCorrectMatrix() {
		double[] eq = new double[mass.length * mass.length];
		for (int i = mass.length; i< eq.length; i++) {
			eq[i] = 1;
		}
		LinearObjectiveFunction f = new LinearObjectiveFunction(eq,0);
		ArrayList<LinearConstraint> constraints = new ArrayList<LinearConstraint>();
		for (int i = 1; i <= mass.length - 1; i++ ) {
			for (int j = i+1; j <= mass.length; j++) {
				double[] eq1 = new double[eq.length];
				double[] eq2 = new double[eq.length];
				eq1[i-1] = 1;
				eq1[j-1] = -1;
				eq2[i-1] = 1;
				eq2[j-1] = -1;
				eq1[mass.length+((2*mass.length-i)*(i-1))/2+(j-i)-1] = 1;				
				eq2[mass.length*(mass.length+1)/2+((2*mass.length-i)*(i-1))/2+(j-i)-1] = -1;
				double[] eq3 = new double[eq.length];
				double[] eq4 = new double[eq.length];
				eq3[mass.length+((2*mass.length-i)*(i-1))/2+(j-i)-1] = 1;				
				eq4[mass.length*(mass.length+1)/2+((2*mass.length-i)*(i-1))/2+(j-i)-1] = 1;
				constraints.add(new LinearConstraint(eq1, Relationship.GEQ, Math.log(mass[i-1][j-1].L)));
		        constraints.add(new LinearConstraint(eq2, Relationship.LEQ, Math.log(mass[i-1][j-1].U)));
		        constraints.add(new LinearConstraint(eq3, Relationship.GEQ, 0));
		        constraints.add(new LinearConstraint(eq4, Relationship.GEQ, 0));
			}
		}
		for (int i = 0; i < mass.length; i++ ) {
			double[] eq1 = new double[eq.length];
			eq1[i] = 1;
			//constraints.add(new LinearConstraint(eq1, Relationship.LEQ, 0));
		}
		
		SimplexSolver solver = new SimplexSolver();
        PointValuePair solution = solver.optimize(f, new LinearConstraintSet(constraints),
                                                  GoalType.MINIMIZE);
        double[] answer = solution.getPoint();
        for (int i = 1; i <= mass.length - 1; i++ ) {
			for (int j = i+1; j <= mass.length; j++) {
				mass[i-1][j-1].L *= Math.exp(answer[mass.length+((2*mass.length-i)*(i-1))/2+(j-i)-1]);
				mass[i-1][j-1].U *= Math.exp(answer[mass.length*(mass.length+1)/2+((2*mass.length-i)*(i-1))/2+(j-i)-1]);
				mass[j-1][i-1].L = 1./mass[i-1][j-1].U;
				mass[j-1][i-1].U = 1./mass[i-1][j-1].L;
			}
		}        
        
	}
	
	public Interval FindWi(int i0) {
		double[] eq = new double[mass.length];
		eq[i0] = 1;
		
		LinearObjectiveFunction f = new LinearObjectiveFunction(eq,0);

		ArrayList<LinearConstraint> constraints = new ArrayList<LinearConstraint>();
		for (int i = 1; i <= mass.length - 1; i++ ) {
			for (int j = i+1; j <= mass.length; j++) {
				double[] eq1 = new double[eq.length];
				double[] eq2 = new double[eq.length];
				eq1[i-1] = 1;
				eq1[j-1] = -mass[i-1][j-1].L;
				eq2[i-1] = 1;
				eq2[j-1] = -mass[i-1][j-1].U;
				constraints.add(new LinearConstraint(eq1, Relationship.GEQ, 0));
		        constraints.add(new LinearConstraint(eq2, Relationship.LEQ, 0));
			}
		}
		
		for (int i = 0; i < mass.length; i++ ) {
			double[] eq1 = new double[eq.length];
			eq1[i] = 1;
			constraints.add(new LinearConstraint(eq1, Relationship.GEQ, 0));
			double[] eq2 = new double[eq.length];
			eq2[i] = 1;
			constraints.add(new LinearConstraint(eq1, Relationship.LEQ, 1));
		}
		double[] eq1 = new double[mass.length];
		for (int i = 0; i < mass.length; i++ ) {
			eq1[i]=1;
		}
		constraints.add(new LinearConstraint(eq1, Relationship.EQ, 1));
		
		SimplexSolver solver = new SimplexSolver();
		double[] ans = solver.optimize(f, new LinearConstraintSet(constraints),
                GoalType.MAXIMIZE).getPoint();
		double sum = 0;
		for (int i = 0; i <ans.length; i++) {
			sum+=ans[i];
		}
		double maxX = ans[i0]/sum;
		ans = solver.optimize(f, new LinearConstraintSet(constraints),
                GoalType.MINIMIZE).getPoint();
		sum = 0;
		for (int i = 0; i <ans.length; i++) {
			sum+=ans[i];
		}
        double minX = ans[i0]/sum;
        
        return new Interval(minX, maxX);      
	}
	
	public Collection<Interval> findW() {
		ArrayList<Interval> w = new ArrayList<>();
		MatrixIsCorrect();
		for (int i = 0; i < mass.length; i++) {
			w.add(FindWi(i));
		}
		System.out.println(name +" w: "+ System.getProperty("line.separator") + w);
		return w;
	}
	
	public String toString() {
		StringBuilder res= new StringBuilder(name + ": "+ System.getProperty("line.separator"));
		for (int i = 0; i < mass.length; i++ ) {
			res.append("[");
			for (int j = 0; j < mass.length; j++) {
				res.append(mass[i][j].toString() + ",");
			}
			res.append("]" + System.getProperty("line.separator"));			
		}
		return res.toString();
	}
}