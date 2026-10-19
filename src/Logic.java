import java.util.Locale;
import java.util.Scanner;

public class Logic {

    private static final String OUTPUT_FILE = "output/output5.txt";

    private int intNumSeg;
    private double dblE;
    private int intDOF;
    private double dblX;
    private double dblInitialX;

    private final SimpsonIntegration simpsonIntegration;
    private final Output output;

    public Logic() {
        simpsonIntegration = new SimpsonIntegration();
        output = new Output();
        dblE = simpsonIntegration.getDblE();
    }

    public void logica() {
        Scanner scanner = new Scanner(System.in);

        try {
            dblInitialX = readDouble(scanner, "Valor inicial: ");
            dblX = readDouble(scanner, "Valor final: ");
            intDOF = readPositiveInt(scanner, "DOF (grados de libertad): ");

            double finalValue = calcular(dblInitialX, dblX, intDOF);
            String result = createResult(finalValue);

            System.out.print(result);
            output.writeData(OUTPUT_FILE, result);
            System.out.println("Resultado guardado en " + OUTPUT_FILE);
        } catch (IllegalArgumentException | IllegalStateException e) {
            System.out.println("No se pudo realizar la integración: " + e.getMessage());
        }
    }

    public double calcular(double initialX, double finalX, int dof) {
        dblInitialX = initialX;
        dblX = finalX;
        intDOF = dof;

        double finalValue = simpsonIntegration.integrate(initialX, finalX, dof);
        intNumSeg = simpsonIntegration.getIntNumSeg();
        dblE = simpsonIntegration.getDblE();
        return finalValue;
    }

    public int getIntNumSeg() {
        return intNumSeg;
    }

    public double getDblE() {
        return dblE;
    }

    private double readDouble(Scanner scanner, String prompt) {
        System.out.print(prompt);
        if (!scanner.hasNextLine()) {
            throw new IllegalStateException("Falta un dato de entrada.");
        }

        String value = scanner.nextLine().trim().replace(',', '.');
        try {
            double parsedValue = Double.parseDouble(value);
            if (!Double.isFinite(parsedValue)) {
                throw new NumberFormatException();
            }
            return parsedValue;
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Se esperaba un número válido.");
        }
    }

    private int readPositiveInt(Scanner scanner, String prompt) {
        System.out.print(prompt);
        if (!scanner.hasNextLine()) {
            throw new IllegalStateException("Falta un dato de entrada.");
        }

        try {
            int value = Integer.parseInt(scanner.nextLine().trim());
            if (value <= 0) {
                throw new NumberFormatException();
            }
            return value;
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("DOF debe ser un entero mayor que cero.");
        }
    }

    private String createResult(double finalValue) {
        return String.format(Locale.US, "%.5f%n", finalValue);
    }
}
