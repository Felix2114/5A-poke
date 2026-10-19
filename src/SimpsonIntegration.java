import java.util.Arrays;

/**
 * Integra la función de densidad t de Student mediante la regla de Simpson.
 */
public class SimpsonIntegration {

    public static final int INITIAL_NUM_SEGMENTS = 10;
    public static final double DEFAULT_ERROR = 0.00001;

    private int intNumSeg;
    private double dblW;
    private double dblE;
    private int intDOF;
    private double dblX;
    private double dblInitialX;
    private double[] dblTotXi;
    private double[] dblFirstBaseTerms;
    private double dblXExponent;
    private double dblCoeff;
    private double[] dblFx;
    private double[] dblFinalTerms;
    private double dblFinalValue;

    private final GammaFunction gammaFunction;

    public SimpsonIntegration() {
        this(DEFAULT_ERROR);
    }

    public SimpsonIntegration(double dblE) {
        if (!Double.isFinite(dblE) || dblE <= 0.0) {
            throw new IllegalArgumentException("El error aceptable debe ser positivo.");
        }

        this.dblE = dblE;
        gammaFunction = new GammaFunction();
    }

    /**
     * Obtiene el área bajo la distribución t entre los dos límites indicados.
     * Se duplica el número de segmentos hasta que dos aproximaciones
     * consecutivas difieren menos que el error aceptable.
     */
    public double integrate(double initialX, double finalX, int intDOF) {
        validateInput(initialX, finalX, intDOF);

        dblInitialX = initialX;
        dblX = finalX;
        this.intDOF = intDOF;

        intNumSeg = INITIAL_NUM_SEGMENTS;
        double previousValue = calculateCurrentApproximation();

        while (true) {
            intNumSeg *= 2;
            double currentValue = calculateCurrentApproximation();

            if (Math.abs(currentValue - previousValue) < dblE) {
                return currentValue;
            }

            previousValue = currentValue;
        }
    }

    /** Calcula W = (valor final - valor inicial) / número de segmentos. */
    public double computeW(int intNumSeg, double dblX) {
        if (intNumSeg <= 0 || intNumSeg % 2 != 0) {
            throw new IllegalArgumentException("El número de segmentos debe ser positivo y par.");
        }

        this.intNumSeg = intNumSeg;
        this.dblX = dblX;
        dblW = (dblX - dblInitialX) / intNumSeg;
        return dblW;
    }

    /** Calcula xi = valor inicial + iW, para i desde 0 hasta numSeg. */
    public double[] computeXi(int intNumSeg) {
        if (intNumSeg != this.intNumSeg) {
            throw new IllegalArgumentException("El número de segmentos no coincide con W.");
        }

        dblTotXi = new double[intNumSeg + 1];
        for (int index = 0; index <= intNumSeg; index++) {
            dblTotXi[index] = dblInitialX + index * dblW;
        }
        return Arrays.copyOf(dblTotXi, dblTotXi.length);
    }

    /** Calcula 1 + (xi² / dof), la base de la potencia de la distribución t. */
    public double[] computeFirstBaseTerms(
        int intNumSeg,
        double[] dblTotXi,
        int intDOF
    ) {
        if (dblTotXi == null || dblTotXi.length != intNumSeg + 1 || intDOF <= 0) {
            throw new IllegalArgumentException("Los datos para calcular la distribución t no son válidos.");
        }

        this.intDOF = intDOF;
        dblFirstBaseTerms = new double[intNumSeg + 1];
        for (int index = 0; index <= intNumSeg; index++) {
            dblFirstBaseTerms[index] = 1.0 + (dblTotXi[index] * dblTotXi[index] / intDOF);
        }
        return Arrays.copyOf(dblFirstBaseTerms, dblFirstBaseTerms.length);
    }

    /** Calcula -(dof + 1) / 2. */
    public double computeXExponent(int intDOF) {
        if (intDOF <= 0) {
            throw new IllegalArgumentException("Los grados de libertad deben ser positivos.");
        }

        this.intDOF = intDOF;
        dblXExponent = -(intDOF + 1.0) / 2.0;
        return dblXExponent;
    }

    /** Calcula el coeficiente Γ((dof + 1)/2) / (sqrt(dofπ) Γ(dof/2)). */
    public double computeCoeff(int intDOF) {
        if (intDOF <= 0) {
            throw new IllegalArgumentException("Los grados de libertad deben ser positivos.");
        }

        this.intDOF = intDOF;
        double numerator = gammaFunction.computeDblGamma((intDOF + 1.0) / 2.0);
        double denominator = Math.sqrt(intDOF * Math.PI)
            * gammaFunction.computeDblGamma(intDOF / 2.0);
        dblCoeff = numerator / denominator;
        return dblCoeff;
    }

    /** Evalúa la función de densidad t en todos los xi. */
    public double[] computeFx(int intNumSeg) {
        if (dblFirstBaseTerms == null || dblFirstBaseTerms.length != intNumSeg + 1) {
            throw new IllegalStateException("Primero deben calcularse los términos base.");
        }

        dblFx = new double[intNumSeg + 1];
        for (int index = 0; index <= intNumSeg; index++) {
            dblFx[index] = dblCoeff * Math.pow(dblFirstBaseTerms[index], dblXExponent);
        }
        return Arrays.copyOf(dblFx, dblFx.length);
    }

    /** Aplica los multiplicadores 1, 4, 2, ..., 4, 1 de Simpson. */
    public double[] computeFinalTerms(int intNumSeg) {
        if (dblFx == null || dblFx.length != intNumSeg + 1) {
            throw new IllegalStateException("Primero deben calcularse los valores de f(xi).");
        }

        dblFinalTerms = new double[intNumSeg + 1];
        for (int index = 0; index <= intNumSeg; index++) {
            if (index == 0 || index == intNumSeg) {
                dblFinalTerms[index] = dblFx[index];
            } else if (index % 2 == 0) {
                dblFinalTerms[index] = 2.0 * dblFx[index];
            } else {
                dblFinalTerms[index] = 4.0 * dblFx[index];
            }
        }
        return Arrays.copyOf(dblFinalTerms, dblFinalTerms.length);
    }

    /** Suma los términos ponderados y los multiplica por W / 3. */
    public double computeFinalValue(int intNumSeg) {
        if (dblFinalTerms == null || dblFinalTerms.length != intNumSeg + 1) {
            throw new IllegalStateException("Primero deben calcularse los términos finales.");
        }

        double sum = 0.0;
        for (double term : dblFinalTerms) {
            sum += term;
        }

        dblFinalValue = dblW * sum / 3.0;
        return dblFinalValue;
    }

    public int getIntNumSeg() {
        return intNumSeg;
    }

    public double getDblW() {
        return dblW;
    }

    public double getDblE() {
        return dblE;
    }

    public double getDblFinalValue() {
        return dblFinalValue;
    }

    private double calculateCurrentApproximation() {
        computeW(intNumSeg, dblX);
        computeXi(intNumSeg);
        computeFirstBaseTerms(intNumSeg, dblTotXi, intDOF);
        computeXExponent(intDOF);
        computeCoeff(intDOF);
        computeFx(intNumSeg);
        computeFinalTerms(intNumSeg);
        return computeFinalValue(intNumSeg);
    }

    private void validateInput(double initialX, double finalX, int intDOF) {
        if (!Double.isFinite(initialX) || !Double.isFinite(finalX)) {
            throw new IllegalArgumentException("Los límites de integración deben ser finitos.");
        }
        if (intDOF <= 0) {
            throw new IllegalArgumentException("Los grados de libertad deben ser mayores que cero.");
        }
    }
}
