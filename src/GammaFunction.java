public class GammaFunction {
    private double gammaValue;

    public double calcularGammaEntero(int intValue) {
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

    public double calcularGammaDecimal(double doubleValue) {
        if (!Double.isFinite(doubleValue) || doubleValue <= 0.0) {
            throw new IllegalArgumentException("El argumento de gamma debe ser positivo y finito.");
        }

        if (doubleValue == Math.rint(doubleValue) && doubleValue <= 171.0) {
            return calcularGammaEntero((int) Math.rint(doubleValue));
        }

        double value = doubleValue;
        double divisor = 1.0;

        while (value < 20.0) {
            divisor = divisor * value;
            value = value + 1.0;
        }

        gammaValue = Math.exp(calcularLogGamma(value)) / divisor;
        return gammaValue;
    }

    public double obtenerValorGamma() {
        return gammaValue;
    }

    private double calcularLogGamma(double value) {
        double correction = 1.0 / (12.0 * value);
        correction = correction - 1.0 / (360.0 * Math.pow(value, 3));
        correction = correction + 1.0 / (1260.0 * Math.pow(value, 5));
        correction = correction - 1.0 / (1680.0 * Math.pow(value, 7));

        return (value - 0.5) * Math.log(value)
            - value
            + 0.5 * Math.log(2.0 * Math.PI)
            + correction;
    }
}
