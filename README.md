# Quantum-Classification

The bellow ExampleCode.mlx is intended to create a Quantum Supervised Classification (QSC) network, while n is the overall number of nodes, k is the number of input nodes (and the size of the Hilbert space of the input states) and q is the number of output nodes (= number of classes to classify to).
Determining this 3 parameters, along with the input states (psi) and the PSO parameters, will output a well trained QSC that can distinguish between different classes.
2 more codes you can use and manage are parallelCF.m and NeuralNetwork.m .
parallelCF.m is the cost function, and you should tune it as you need and as described inside the script.
NeuralNetwork.m is the basic code of the QSC network, working on the basis of the Lindblad equation. you can use it to test the trained network, or to see how iy works.
For any other questions please contact the corresponding author.

Enjoy.
