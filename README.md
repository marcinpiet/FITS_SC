# FITS_SC: SuperCollider Classes for FITS File Handling

This repository provides SuperCollider classes designed for comprehensive interaction with FITS (Flexible Image Transport System) files. FITS is a standard digital file format for astronomical data, and these classes aim to bring its capabilities directly into the SuperCollider environment for sound synthesis, analysis, and visualization.

## Features

*   **FITSFile.sc:**
    *   **Reading:** Load FITS files and access their header information (keywords, values, comments).
    *   **Data Access:** Retrieve image data (primary data unit or extensions) as SuperCollider arrays or buffers.
    *   **Metadata:** Query and manipulate FITS header metadata.
    *   **Editing (Planned/Future):** Capabilities for modifying FITS header keywords and data.

*   **FITSVisualisation.sc:**
    *   **Data Mapping:** Map FITS image data to sound parameters (e.g., pixel intensity to amplitude, frequency, or filter cutoff).
    *   **Real-time Sonification:** Tools for sonifying FITS data in real-time within SuperCollider.
    *   **Visual Representation (Planned/Future):** Integration with SuperCollider's GUI or external visualization tools to display FITS images.

## Getting Started

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/marcinpiet/FITS_SC.git
    ```
2.  **Place in SuperCollider Extensions:**
    Create a symbolic link from the cloned directory to your SuperCollider Extensions folder:
    ```bash
    ln -s /path/to/your/cloned/FITS_SC /Users/yourusername/Library/Application\ Support/SuperCollider/Extensions/FITS_SC
    ```
    (Replace `/path/to/your/cloned/FITS_SC` with the actual path where you cloned this repository.)
3.  **Recompile Class Library:** In SuperCollider, go to `Language -> Recompile Class Library`.

## Contribution

Contributions are welcome! Please feel free to open issues or pull requests.

## License

[Add License Information Here, e.g., MIT, GPL, etc.]