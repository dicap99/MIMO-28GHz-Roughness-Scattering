# MIMO-28GHz-Roughness-Scattering

Este repositorio contiene los datos experimentales y los algoritmos de simulación en MATLAB desarrollados para el estudio de la viabilidad de multiplexación espacial en sistemas MIMO operando en la banda milimétrica (28 GHz) sobre superficies rugosas.

## Descripción del Proyecto
El objetivo principal de esta investigación fue evaluar cómo la rugosidad superficial (específicamente lechos de grava y superficies dieléctricas) afecta el comportamiento del canal de dos trayectorias a nivel electromagnético y polarimétrico. 

Mediante mediciones empíricas y modelado matemático (Rayleigh-Ament e IEM/GO calibrado), se demostró cómo la dispersión difusa genera una pérdida drástica en el aislamiento de polarización cruzada (XPD). A partir de la matriz del canal reconstruida empíricamente, se implementa la Descomposición de Valores Singulares (SVD) y el algoritmo de **Waterfilling** para cuantificar los Grados de Libertad (DoF). El estudio concluye demostrando analíticamente cómo un canal MIMO $4 \times 4$ colapsa a un esquema asimétrico equivalente a SISO ante la alta rugosidad.

## Estructura del Repositorio

- **`/Mediciones`**: Carpeta que contiene la base de datos completa con los archivos en formato CSV resultantes de la campaña experimental. Los datos incluyen las mediciones polarimétricas ($P_{vv}, P_{vh}, P_{hv}, P_{hh}$) bajo diversos ángulos de incidencia geométrica ($30^\circ$, $45^\circ$, $60^\circ$) frente a escenarios de referencia lisa, grava fina y grava gruesa.
- **`Waterfilling_30_45_60.m`**: Script principal de MATLAB que extrae la matriz de canal empírica a partir de los datos recolectados, aplica descomposición SVD y ejecuta el algoritmo de Waterfilling para visualizar gráficamente el colapso de la capacidad en los escenarios evaluados.
- **`Waterfilling_Empirico_Simulacion.m`**: Script de procesamiento avanzado y validación matemática enfocado en simular las curvas continuas de capacidad frente a rugosidad.

## Requisitos
- MATLAB (R2021a o superior recomendado).
- Toolbox básicos de análisis matricial.

## Instrucciones de Uso
1. Clona o descarga este repositorio en tu ordenador local.
2. Abre MATLAB y navega hacia la carpeta raíz del repositorio como tu entorno de trabajo (*Current Folder*).
3. Asegúrate de que los archivos de la carpeta `/Mediciones` mantengan su ruta relativa para que los algoritmos puedan leer los `.csv` correctamente.
4. Ejecuta el script `Waterfilling_30_45_60.m` para generar de manera automática las gráficas de potencia dispersada, capacidad del canal y distribución óptima de energía en los subcanales.
