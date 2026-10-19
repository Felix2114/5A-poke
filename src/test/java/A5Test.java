import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.InputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;


public class A5Test {

    private static final double TOLERANCE = 0.00001;

    @Test
    void calculaGammaParaValoresEnterosYNoEnteros() {
        GammaFunction gamma = new GammaFunction();

        assertEquals(1.0, gamma.calcularGammaEntero(1));
        assertEquals(24.0, gamma.calcularGammaEntero(5));
        assertEquals(Math.sqrt(Math.PI), gamma.calcularGammaDecimal(0.5), 1e-12);
        assertEquals(1.329340388179137, gamma.calcularGammaDecimal(2.5), 1e-12);
        assertEquals(3.625609908221908, gamma.calcularGammaDecimal(0.25), 1e-12);
        assertEquals(3.625609908221908, gamma.obtenerValorGamma(), 1e-12);
        assertEquals(
            Double.POSITIVE_INFINITY,
            gamma.calcularGammaDecimal((double) Integer.MAX_VALUE + 1.0)
        );
    }

    @Test
    void rechazaArgumentosInvalidosDeGamma() {
        GammaFunction gamma = new GammaFunction();

        assertThrows(IllegalArgumentException.class, () -> gamma.calcularGammaEntero(0));
        assertThrows(IllegalArgumentException.class, () -> gamma.calcularGammaDecimal(0.0));
        assertThrows(IllegalArgumentException.class, () -> gamma.calcularGammaDecimal(Double.NaN));
    }

    @Test
    void calculaLosTresCasosDeReferenciaDelPdf() {
        SimpsonIntegration integration = new SimpsonIntegration();

        assertEquals(0.30006, integration.integrate(0.0, 1.1, 9), TOLERANCE);
        assertEquals(0.36757, integration.integrate(0.0, 1.1812, 10), TOLERANCE);
        assertEquals(0.49500, integration.integrate(0.0, 2.750, 30), TOLERANCE);
    }

    @Test
    void refinaLosSegmentosCuandoSeSolicitaMayorPrecision() {
        SimpsonIntegration integration = new SimpsonIntegration(1e-12);

        assertEquals(0.3500586190716141, integration.integrate(0.0, 1.1, 9), 1e-12);
        assertTrue(integration.getIntNumSeg() > 20);
    }

    @Test
    void calculaLosTerminosDeLaReglaDeSimpson() {
        SimpsonIntegration integration = new SimpsonIntegration();

        assertEquals(0.1, integration.computeW(10, 1.0), 0.0);
        assertEquals(0.1, integration.getDblW(), 0.0);
        assertArrayEquals(
            new double[] {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0},
            integration.computeXi(10),
            1e-12
        );

        double[] bases = integration.computeFirstBaseTerms(
            10,
            integration.computeXi(10),
            9
        );
        assertEquals(1.0, bases[0], 0.0);
        assertEquals(1.0 + (1.0 / 9.0), bases[10], 1e-12);
        assertEquals(-5.0, integration.computeXExponent(9), 0.0);
        assertEquals(0.3880349, integration.computeCoeff(9), 0.000001);

        double[] functionValues = integration.computeFx(10);
        assertEquals(0.3880349, functionValues[0], 0.000001);

        double[] finalTerms = integration.computeFinalTerms(10);
        assertEquals(functionValues[0], finalTerms[0], 0.0);
        assertEquals(4.0 * functionValues[1], finalTerms[1], 0.0);
        assertEquals(2.0 * functionValues[2], finalTerms[2], 0.0);
        assertEquals(functionValues[10], finalTerms[10], 0.0);
        assertTrue(integration.computeFinalValue(10) > 0.0);
        assertEquals(integration.computeFinalValue(10), integration.getDblFinalValue(), 0.0);
    }

    @Test
    void validaLosDatosYElOrdenDeLosCalculos() {
        SimpsonIntegration integration = new SimpsonIntegration();

        assertThrows(IllegalArgumentException.class, () -> new SimpsonIntegration(0.0));
        assertThrows(IllegalArgumentException.class, () -> new SimpsonIntegration(Double.NaN));
        assertThrows(IllegalArgumentException.class, () -> integration.integrate(Double.NaN, 1.0, 9));
        assertThrows(IllegalArgumentException.class, () -> integration.integrate(0.0, Double.NaN, 9));
        assertThrows(IllegalArgumentException.class, () -> integration.integrate(0.0, 1.0, 0));
        assertThrows(IllegalArgumentException.class, () -> integration.computeW(0, 1.0));
        assertThrows(IllegalArgumentException.class, () -> integration.computeW(3, 1.0));

        integration.computeW(10, 1.0);
        assertThrows(IllegalArgumentException.class, () -> integration.computeXi(2));
        assertThrows(
            IllegalArgumentException.class,
            () -> integration.computeFirstBaseTerms(10, null, 9)
        );
        assertThrows(
            IllegalArgumentException.class,
            () -> integration.computeFirstBaseTerms(10, new double[2], 9)
        );
        assertThrows(
            IllegalArgumentException.class,
            () -> integration.computeFirstBaseTerms(10, new double[11], 0)
        );
        assertThrows(IllegalArgumentException.class, () -> integration.computeXExponent(0));
        assertThrows(IllegalArgumentException.class, () -> integration.computeCoeff(0));

        SimpsonIntegration noCalculatedTerms = new SimpsonIntegration();
        assertThrows(IllegalStateException.class, () -> noCalculatedTerms.computeFx(10));
        assertThrows(IllegalStateException.class, () -> noCalculatedTerms.computeFinalTerms(10));
        assertThrows(IllegalStateException.class, () -> noCalculatedTerms.computeFinalValue(10));

        integration.computeFirstBaseTerms(10, integration.computeXi(10), 9);
        integration.computeXExponent(9);
        integration.computeCoeff(9);
        integration.computeFx(10);
        assertThrows(IllegalStateException.class, () -> integration.computeFx(2));
        assertThrows(IllegalStateException.class, () -> integration.computeFinalTerms(2));
        integration.computeFinalTerms(10);
        assertThrows(IllegalStateException.class, () -> integration.computeFinalValue(2));
    }

    @Test
    void logicCalculaEnUnIntervaloConValorInicialDistintoDeCero() {
        Logic logic = new Logic();

        assertEquals(0.70011728, logic.calcular(-1.1, 1.1, 9), TOLERANCE);
        assertTrue(logic.getIntNumSeg() >= SimpsonIntegration.INITIAL_NUM_SEGMENTS);
        assertEquals(SimpsonIntegration.DEFAULT_ERROR, logic.getDblE(), 0.0);
    }

    @Test
    void appSolicitaLosTresDatosYElArchivoContieneSoloElResultado() throws Exception {
        new App();
        String console = executeWithInput("0\n1,1\n9\n", () -> App.main(new String[0]));

        assertTrue(console.contains("Valor inicial:"));
        assertTrue(console.contains("Valor final:"));
        assertTrue(console.contains("DOF (grados de libertad):"));
        assertEquals("0.35006\n", Files.readString(Path.of("output", "output5.txt")));
    }

    @Test
    void logicMuestraUnErrorParaUnDofInvalido() {
        String console = executeWithInput("0\n1.1\n0\n", () -> new Logic().logica());

        assertTrue(console.contains("DOF debe ser un entero mayor que cero."));
    }

    @Test
    void logicManejaEntradasFaltantesYValoresNoValidos() {
        assertTrue(
            executeWithInput("", () -> new Logic().logica())
                .contains("Falta un dato de entrada.")
        );
        assertTrue(
            executeWithInput("texto\n", () -> new Logic().logica())
                .contains("Se esperaba un número válido.")
        );
        assertTrue(
            executeWithInput("NaN\n", () -> new Logic().logica())
                .contains("Se esperaba un número válido.")
        );
        assertTrue(
            executeWithInput("0\n1.1\n", () -> new Logic().logica())
                .contains("Falta un dato de entrada.")
        );
    }

    @Test
    void outputEscribeEnUnDirectorioQueNoExiste() throws Exception {
        Path directory = Files.createTempDirectory("simpson-output-");
        Path file = directory.resolve("subdirectorio").resolve("resultado.txt");

        new Output().writeData(file.toString(), "0.35006\n");

        assertEquals("0.35006\n", Files.readString(file));
        Files.deleteIfExists(file);
        Files.deleteIfExists(file.getParent());
        Files.deleteIfExists(directory);
    }

    @Test
    void outputManejaRutasSinDirectorioYErroresDeEscritura() throws Exception {
        Path rootFile = Path.of("resultado-raiz-prueba.txt");
        Output output = new Output();

        output.writeData(rootFile.toString(), "0.49500\n");
        assertEquals("0.49500\n", Files.readString(rootFile));
        Files.deleteIfExists(rootFile);

        Path directory = Files.createTempDirectory("simpson-output-error-");
        String console = executeWithInput(
            "",
            () -> output.writeData(directory.toString(), "no se puede escribir sobre un directorio")
        );
        assertTrue(console.contains("Error al escribir"));
        Files.deleteIfExists(directory);
    }

    private String executeWithInput(String input, Runnable action) {
        InputStream originalInput = System.in;
        PrintStream originalOutput = System.out;
        ByteArrayOutputStream capturedOutput = new ByteArrayOutputStream();

        try {
            System.setIn(new ByteArrayInputStream(input.getBytes(StandardCharsets.UTF_8)));
            System.setOut(new PrintStream(capturedOutput, true, StandardCharsets.UTF_8));
            action.run();
            return capturedOutput.toString(StandardCharsets.UTF_8);
        } finally {
            System.setIn(originalInput);
            System.setOut(originalOutput);
        }
    }
}
