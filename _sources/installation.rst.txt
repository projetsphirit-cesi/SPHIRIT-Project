Installation
============

Prérequis
---------

- Windows 10/11
- `Python 3.11+ <https://www.python.org/downloads/>`_
- `Visual Studio Code <https://code.visualstudio.com/>`_
- `ParaView <https://www.paraview.org/download/>`_ (pour la visualisation des résultats)

.. note::
   Lors de l'installation de Python, cochez bien **"Add Python to PATH"** en bas de la première fenêtre.

Extensions VS Code
------------------

Installez les extensions suivantes depuis l'onglet Extensions (``Ctrl+Shift+X``) :

- **Python** (Microsoft)
- **Jupyter** (Microsoft)

Création de l'environnement virtuel
------------------------------------

Ouvrez un terminal **Git Bash** dans VS Code et exécutez :

.. code-block:: bash

   python -m venv C:/venvs/sph
   source C:/venvs/sph/Scripts/activate
   pip install numpy scipy matplotlib

.. warning::
   Ne créez pas le venv dans un dossier contenant des espaces (ex: OneDrive).
   Utilisez un chemin simple comme ``C:/venvs/sph``.

Configuration du kernel Jupyter
---------------------------------

1. Ouvrez le notebook ``.ipynb`` dans VS Code
2. Cliquez sur **"Select Kernel"** en haut à droite
3. Choisissez **"Python Environments"**
4. Sélectionnez ``C:\venvs\sph\Scripts\python.exe``

Vérification
------------

Exécutez la première cellule du notebook. Si aucune erreur n'apparaît, l'installation est réussie.