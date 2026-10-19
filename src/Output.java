import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;


public class Output {

    public void writeData(String outFile, String outText) {
        try {
            Path outputPath = Path.of(outFile);
            Path parentDirectory = outputPath.getParent();
            if (parentDirectory != null) {
                Files.createDirectories(parentDirectory);
            }
            Files.writeString(outputPath, outText);
        } catch (IOException e) {
            System.out.println("Error al escribir " + outFile + ": " + e.getMessage());
        }
    }
}
