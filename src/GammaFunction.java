public class GammaFunction {

    private static final double[] LANCZOS_COEFFICIENTS = {
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7
    };

    private double gammaValue;

    public double computeIntGamma(int intValue) {
        if (intValue <= 0) {
            throw new IllegalArgumentException("Gamma solo está definida aquí para valores positivos.");
        }

        double result = 1.0;
        for (int value = 2; value < intValue; value++) {
            result *= value;
        }

        gammaValue = result;
        return gammaValue;
    }

    public double computeDblGamma(double doubleValue) {
        if (!Double.isFinite(doubleValue) || doubleValue <= 0.0) {
            throw new IllegalArgumentException("El argumento de gamma debe ser positivo y finito.");
        }

        if (doubleValue == Math.rint(doubleValue)
            && doubleValue <= Integer.MAX_VALUE) {
            return computeIntGamma((int) Math.rint(doubleValue));
        }

        gammaValue = Math.exp(logGamma(doubleValue));
        return gammaValue;
    }

    public double getGammaValue() {
        return gammaValue;
    }

    private double logGamma(double value) {
        if (value < 0.5) {
            return Math.log(Math.PI)
                - Math.log(Math.sin(Math.PI * value))
                - logGamma(1.0 - value);
        }

        double adjustedValue = value - 1.0;
        double series = LANCZOS_COEFFICIENTS[0];

        for (int index = 1; index < LANCZOS_COEFFICIENTS.length; index++) {
            series += LANCZOS_COEFFICIENTS[index] / (adjustedValue + index);
        }

        double t = adjustedValue + LANCZOS_COEFFICIENTS.length - 1.5;
        return 0.5 * Math.log(2.0 * Math.PI)
            + (adjustedValue + 0.5) * Math.log(t)
            - t
            + Math.log(series);
    }
}
