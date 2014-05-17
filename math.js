/**
 * Mathematics library for Javascript
 * Version PROTOTYPE
 *
 * Written by Heaven(Cichol Gricenchos)
 * License under the MIT License
 *
 * heaven@isdev.kr
 * http://www.heavenlab.kr/
 **/

(function() {
	"use strict";

	var Mathematics = {
		// Default arithmetic caculation
		Add: function(A, B) {
			A = Mathematics.isComplex(B) && !Mathematics.isComplex(A) ? Mathematics.Complex(A, 0) : A;
			B = Mathematics.isComplex(A) && !Mathematics.isComplex(B) ? Mathematics.Complex(B, 0) : B;

			return Mathematics.isComplex(A) ? A.Add(B) : A + B;
		},
		Subtract: function(A, B) {
			A = Mathematics.isComplex(B) && !Mathematics.isComplex(A) ? Mathematics.Complex(A, 0) : A;
			B = Mathematics.isComplex(A) && !Mathematics.isComplex(B) ? Mathematics.Complex(B, 0) : B;

			return Mathematics.isComplex(A) ? A.Subtract(B) : A - B;
		},
		Multiply: function(A, B) {
			A = Mathematics.isComplex(B) && !Mathematics.isComplex(A) ? Mathematics.Complex(A, 0) : A;
			B = Mathematics.isComplex(A) && !Mathematics.isComplex(B) ? Mathematics.Complex(B, 0) : B;

			return Mathematics.isComplex(A) ? A.Multiply(B) : A * B;
		},
		Divide: function(A, B) {
			A = Mathematics.isComplex(B) && !Mathematics.isComplex(A) ? Mathematics.Complex(A, 0) : A;
			B = Mathematics.isComplex(A) && !Mathematics.isComplex(B) ? Mathematics.Complex(B, 0) : B;

			return Mathematics.isComplex(A) ? A.Divide(B) : A / B;
		},

		// Essential algebric functions
		Power: function(W, Z) {
			return Mathematics.isComplex(W) || Mathematics.isComplex(Z) ? Mathematics.Exp(Mathematics.Multiply(Z, Mathematics.Ln(W))) : Math.pow(W, Z);
		},
		Root: function(X, N) {
			if(!N) {
				N = 2;
			}

			return Mathematics.Power(X, Mathematics.Divide(1, N));
		},
		Square: function(X) {
			return Mathematics.Power(X, 2);
		},
		Cube: function(X) {
			return Mathematics.Power(X, 3);
		},
		Abs: function(X) {
			return Mathematics.isComplex(X) ? Mathematics.Distance([X.RealPart(), X.ImaginaryPart()], [0, 0]) : Math.abs(X);
		},
		Sign: function(X) {
			return X == 0 ? X : X / Mathematics.Abs(X);
		},
		Floor: function(X) {
			return Mathematics.isComplex(X) ? Mathematics.Complex(Mathematics.Floor(X.RealPart()), Mathematics.Floor(X.ImaginaryPart())) : Math.floor(X);
		},
		Ceil: function(X) {
			return Mathematics.isComplex(X) ? Mathematics.Complex(Mathematics.Ceil(X.RealPart()), Mathematics.Ceil(X.ImaginaryPart())) : Math.ceil(X);	
		},
		Round: function(X) {
			return Mathematics.isComplex(X) ? Mathematics.Complex(Mathematics.Round(X.RealPart()), Mathematics.Round(X.ImaginaryPart())) : Math.round(X);	
		},
		GCD: function(A, B) {
			A = Mathematics.Abs(A);
			B = Mathematics.Abs(B);

			if(B > A) {
				return Mathematics.GCD(B, A);
			}

			if(A % B == 0) {
				return B;
			}

			return Mathematics.GCD(B, A % B);
		},
		LCM: function(A, B) {
			return A * B / Mathematics.GCD(A, B);
		},

		// Trigonal functions
		Sin: function(X) {
			return Mathematics.isComplex(X) ? Mathematics.Complex(Mathematics.Multiply(Mathematics.Sin(X.RealPart()), Mathematics.Cosh(X.ImaginaryPart())), Mathematics.Multiply(Mathematics.Cos(X.RealPart()), Mathematics.Sinh(X.ImaginaryPart()))) : Math.sin(X);
		},
		Cos: function(X) {
			return Mathematics.isComplex(X) ? Mathematics.Complex(Mathematics.Multiply(Mathematics.Cos(X.RealPart()), Mathematics.Cosh(X.ImaginaryPart())), Mathematics.Multiply(Mathematics.Sin(X.RealPart()), Mathematics.Sinh(X.ImaginaryPart()))).Conjugate() : Math.cos(X);
		},
		Tan: function(X) {
			return Mathematics.Divide(Mathematics.Sin(X), Mathematics.Cos(X));
		},
		Sec: function(X) {
			return Mathematics.Divide(1, Mathematics.Cos(X));
		},
		Csc: function(X) {
			return Mathematics.Divide(1, Mathematics.Sin(X));
		},
		Cot: function(X) {
			return Mathematics.Divide(1, Mathematics.Tan(X));
		},
		Arcsin: function(X) {
			X = Mathematics.isComplex(X) ? X : Mathematics.Complex(X, 0);
			
			return Mathematics.Complex(0, 1).Multiply(Mathematics.Arcsinh(Mathematics.Complex(X.ImaginaryPart(), X.RealPart())));
		},
		Arccos: function(X) {
			return Mathematics.Subtract(Mathematics.PI / 2, Mathematics.Arcsin(X));
		},
		Arctan: function(X) {
			return Mathematics.Arcsin(Mathematics.Divide(X, Mathematics.Root(Mathematics.Add(Mathematics.Square(X), 1))));
		},
		Arcsec: function(X) {
			return Mathematics.Arccos(Mathematics.Divide(1, X));
		},
		Arccsc: function(X) {
			return Mathematics.Arcsin(Mathematics.Divide(1, X));
		},
		Arccot: function(X) {
			return Mathematics.Arctan(Mathematics.Divide(1, X));
		},

		// Hyperbolic functions
		Sinh: function(X) {
			return Mathematics.Divide(Mathematics.Subtract(Mathematics.Exp(X), Mathematics.Exp(Mathematics.Multiply(X, -1))), 2);
		},
		Cosh: function(X) {
			return Mathematics.Divide(Mathematics.Add(Mathematics.Exp(X), Mathematics.Exp(Mathematics.Multiply(X, -1))), 2);
		},
		Tanh: function(X) {
			return Mathematics.Divide(Mathematics.Sinh(X), Mathematics.Cosh(X));
		},
		Sech: function(X) {
			return Mathematics.Divide(1, Mathematics.Cosh(X));
		},
		Csch: function(X) {
			return Mathematics.Divide(1, Mathematics.Sinh(X));	
		},
		Coth: function(X) {
			return Mathematics.Divide(1, Mathematics.Tanh(X));
		},
		Arcsinh: function(X) {
			return Mathematics.Ln(Mathematics.Add(X, Mathematics.Root(Mathematics.Add(Mathematics.Square(X), 1))));
		},
		Arccosh: function(X) {
			return Mathematics.Ln(Mathematics.Add(X, Mathematics.Root(Mathematics.Subtract(Mathematics.Square(X), 1))));
		},
		Arctanh: function(X) {
			return Mathematics.Divide(Mathematics.Ln(Mathematics.Divide(Mathematics.Add(1, X), Mathematics.Subtract(1, X))), 2);
		},
		Arcsech: function(X) {
			return Mathematics.Arccosh(Mathematics.Divide(1, X));
		},
		Arccsch: function(X) {
			return Mathematics.Arcsinh(Mathematics.Divide(1, X));
		},
		Arccoth: function(X) {
			return Mathematics.Arctanh(Mathematics.Divide(1, X));
		},

		// Exponent and logarithm
		Exp: function(X) {
			return Mathematics.isComplex(X) ? Mathematics.Complex(Mathematics.Exp(X.RealPart()) * Mathematics.Cos(X.ImaginaryPart()), Mathematics.Exp(X.RealPart()) * Mathematics.Sin(X.ImaginaryPart())) : Math.exp(X);
		},
		Log: function(X) {
			return Mathematics.isComplex(X) ? Mathematics.Ln(X).Divide(Mathematics.Complex(Mathematics.Ln(10), 0)) : Math.log(X) / Math.log(10);
		},
		Lb: function(X) {
			return Mathematics.isComplex(X) ? Mathematics.Ln(X).Divide(Mathematics.Complex(Mathematics.Ln(2), 0)) : Math.log(X) / Math.log(2);
		},
		Ln: function(X) {
			return Mathematics.isComplex(X) ? Mathematics.Complex(Mathematics.Ln(Mathematics.Abs(X)), Math.atan2(X.ImaginaryPart(), X.RealPart())) : (X < 0 ? Mathematics.Complex(Mathematics.Ln(-X), Mathematics.PI) : Math.log(X));
		},

		// Sequence (Progression)
		AP: function(A, D) {
			return function(N) {
				return A + (N - 1) * D;
			};
		},
		GP: function(A, R) {
			return function(N) {
				return A * Mathematics.Power(A, R - 1);
			};
		},
		HP: function(A, D) {
			return function(N) {
				return 1 / Mathematics.AP(A, D)(N);
			};
		},
		Fibonacci: function(N) {
			return Mathematics.Round((1 / Mathematics.Root(5)) * (Mathematics.Power((1 + Mathematics.Root(5)) / 2, N) - Mathematics.Power((1 - Mathematics.Root(5)) / 2, N)));
		},
		Sum: function(From, To, Fx) {
			var Summation = Mathematics.Complex(0, 0);

			for(var Iterator = From; Iterator <= To; Iterator++) {
				if(Mathematics.Abs(To) == Infinity) {
					if(Summation.Stringify() != Mathematics.Add(Summation, Fx(Iterator)).Stringify()) {
						Summation = Mathematics.Add(Summation, Fx(Iterator));
					} else {
						break;
					}
				} else {
					Summation = Mathematics.Add(Summation, Fx(Iterator));
				}
			}

			return Summation.Realize ? Summation.Realize() : Summation;
		},
		Product: function(From, To, Fx) {
			var Production = Mathematics.Complex(1, 0);

			for(var Iterator = From; Iterator <= To; Iterator++) {
				if(Mathematics.Abs(To) == Infinity) {
					if(Production.Stringify() != Mathematics.Multiply(Production, Fx(Iterator)).Stringify()) {
						Production = Mathematics.Multiply(Production, Fx(Iterator));
					} else {
						break;
					}
				} else {
					Production = Mathematics.Multiply(Production, Fx(Iterator));
				}
			}

			return Production.Realize ? Production.Realize() : Production;
		},

		// Extra functions
		Sinc: function(X) {
			X = Mathematics.isComplex(X) && X.Realize ? X.Realize() : X;

			return X == 0 ? 1 : Mathematics.Divide(Mathematics.Sin(X), X);
		},
		Zeta: function(S) {
			return S > 1 ? (S == Infinity ? 1 : Mathematics.Sum(1, Infinity, function(N) { return Mathematics.Divide(1, Mathematics.Power(N, S)) })) : undefined;
		},
		Gamma: function(T) {
			return Mathematics.isInteger(T) ? Mathematics.Product(1, Mathematics.Subtract(T, 1), function(N) { return N }) : Mathematics.Divide(Mathematics.Product(1, Infinity, function(N) { return Mathematics.Divide(Mathematics.Power(Mathematics.Add(1, Mathematics.Divide(1, N)), T), Mathematics.Add(1, Mathematics.Divide(T, N))) }), T);
		},
		Factorial: function(N) {
			return Mathematics.Gamma(Mathematics.Add(N, 1));
		},
		isInteger: function(X) {
			return Mathematics.isComplex(X) ? (X.Realize instanceof Function && X.Realize() % 1 == 0) : (X % 1 == 0);
		},
		isComplex: function(Z) {
			return (Z.RealPart instanceof Function) && (Z.ImaginaryPart instanceof Function);
		},
		isCoprime: function(A, B) {
			A = Mathematics.isInteger(A) && Mathematics.isComplex(A) ? A.Realize() : A;
			B = Mathematics.isInteger(B) && Mathematics.isComplex(B) ? B.Realize() : B;

			return Mathematics.isInteger(A) && Mathematics.isInteger(B) && Mathematics.GCD(A, B) == 1;
		},

		// Unit conversion
		Radian: function(Degree) {
			return Degree / 180 * Mathematics.PI;
		},
		Degree: function(Radian) {
			return Radian * 180 / Mathematics.PI;
		},

		// Point
		Distance: function(P, Q) {
			return Mathematics.Root(Mathematics.Square(Mathematics.Subtract(Q[0], P[0])) + Mathematics.Square(Mathematics.Subtract(Q[1], P[1])));
		},
		Midpoint: function(P, Q) {
			return [Mathematics.Divide(Mathematics.Add(Q[0], P[0]), 2), Mathematics.Divide(Mathematics.Add(Q[1], P[1]), 2)];
		},

		// Matrix
		Matrix: function(Values) {
			var Matrix = function(Values) {
				var ValueSet = Values, Rows = Values.length, Columns = Values[0].length;

				this.Rows = Rows;
				this.Columns = Columns;

				this.Item = function(Row, Column) {
					if(Row > 0 && Column > 0) {
						return ValueSet[Row - 1][Column - 1];
					}
				};

				this.RowVector = function(Row) {
					var Vector = [];

					for(var ColumnIterator = 1; ColumnIterator <= Columns; ColumnIterator++) {
						Vector.push(this.Item(Row, ColumnIterator));
					}

					return Mathematics.Matrix(Vector);
				};

				this.ColumnVector = function(Column) {
					var Vector = [];

					for(var RowIterator = 1; RowIterator <= Rows; RowIterator++) {
						Vector.push([this.Item(RowIterator, Column)]);
					}

					return Mathematics.Matrix(Vector);
				};

				this.ScalarMultiply = function(Scalar) {
					for(var RowIterator = 0; RowIterator < Rows; RowIterator++) {
						for(var ColumnIterator = 0; ColumnIterator < Columns; ColumnIterator++) {
							ValueSet[RowIterator][ColumnIterator] = Mathematics.Multiply(ValueSet[RowIterator][ColumnIterator], Scalar);
						}
					}
				};

				this.Transpose = function() {
					var Values = [];

					for(var RowIterator = 0; RowIterator < Columns; RowIterator++) {
						Values[RowIterator] = [];

						for(var ColumnIterator = 0; ColumnIterator < Rows; ColumnIterator++) {
							Values[RowIterator][ColumnIterator] = ValueSet[ColumnIterator][RowIterator];
						}
					}

					return Mathematics.Matrix(Values);
				};

				this.Multiply = function(Matrix) {
					if(Columns == Matrix.Rows) {
						var Induction = [], Transposal = Matrix.Transpose();

						for(var RowIterator = 0; RowIterator < Rows; RowIterator++) {
							Induction[RowIterator] = [];

							for(var ColumnIterator = 0; ColumnIterator < Matrix.Columns; ColumnIterator++) {
								Induction[RowIterator][ColumnIterator] = 0;

								for(var Iterator = 1; Iterator <= Columns; Iterator++) {
									Induction[RowIterator][ColumnIterator] = Mathematics.Add(Mathematics.Multiply(this.Item(RowIterator + 1, Iterator), Transposal.Item(ColumnIterator + 1, Iterator)));
								}
							}
						}

						return Mathematics.Matrix(Induction);
					}
				};

				if(Columns == Rows) {
					this.Adjugate = function() {
						var Matrix = [];

						for(var RowIterator = 1; RowIterator <= Rows; RowIterator++) {
							Matrix.push([]);

							for(var ColumnIterator = 1; ColumnIterator <= Columns; ColumnIterator++) {
								Matrix[Matrix.length - 1].push(this.Cofactor(RowIterator, ColumnIterator));
							}
						}

						return Mathematics.Matrix(Matrix).Transpose();
					};

					this.Minor = function(Row, Column) {
						var Minor = [];

						for(var RowIterator = 0; RowIterator < Rows; RowIterator++) {
							if(RowIterator == Row - 1) {
								continue;
							}

							Minor.push([]);

							for(var ColumnIterator = 1; ColumnIterator <= Columns; ColumnIterator++) {
								if(ColumnIterator == Column) {
									continue;
								}

								Minor[Minor.length - 1].push(this.Item(RowIterator + 1, ColumnIterator));
							}
						}

						return Mathematics.Matrix(Minor).Determinant();
					};

					this.Cofactor = function(Row, Column) {
						return Mathematics.Multiply(Mathematics.Power(-1, Row + Column), this.Minor(Row, Column));
					};

					this.Trace = function(Matrix) {
						var Trace = 0;

						for(var Iterator = 1; Iterator <= Rows; Iterator++) {
							Trace = Mathematics.Add(Trace, this.Item(Iterator, Iterator));
						}

						return Trace;
					};

					this.Determinant = function() {
						if(Rows == 1) {
							return this.Item(1, 1);
						} else if(Rows == 2) {
							return this.Item(1, 1) * this.Item(2, 2) - this.Item(1, 2) * this.Item(2, 1);
						} else {
							return Mathematics.Sum(1, Rows, (function(Matrix) {
								return function(I) { return Mathematics.Multiply(Matrix.Item(I, 1), Matrix.Cofactor(I, 1)) }
							})(this));
						}
					};

					if(this.Determinant() != 0) {
						this.InverseMatrix = function() {
							var InverseMatrix = this.Adjugate();
							InverseMatrix.ScalarMultiply(Mathematics.Divide(1, this.Determinant()));

							return InverseMatrix;
						};
					}
				}

				if(Columns == 1 || Rows == 1) {
					this.Serialize = function() {
						if(Columns == 1) {
							return this.Transpose().Serialize();
						} else {
							var SerialArray = [];

							for(var ColumnsIterator = 1; ColumnsIterator <= Columns; ColumnsIterator++) {
								SerialArray.push(this.Item(1, ColumnsIterator));
							}

							return SerialArray;
						}
					};
				}

				this.Stringify = function() {
					var MatrixString = '[';

					for(var RowIterator = 1; RowIterator <= Rows; RowIterator++) {
						MatrixString += RowIterator == 1 ? '' : ' ';

						for(var ColumnIterator = 1; ColumnIterator <= Columns; ColumnIterator++) {
							MatrixString += this.Item(RowIterator, ColumnIterator).Stringify() + (RowIterator == Rows && ColumnIterator == Columns ? '' : ', ');
						}

						MatrixString += RowIterator == Rows ? ']' : '\n';
					}

					return MatrixString;
				};
			};

			return new Matrix(Values);
		},

		// Complex numbers
		Complex: function(A, B) {
			var Complex = function(A, B) {
				var RealPart = A, ImaginaryPart = B;

				this.RealPart = function() {
					return RealPart;
				};

				this.ImaginaryPart = function() {
					return ImaginaryPart;
				};

				if(ImaginaryPart == 0) {
					this.Realize = function() {
						return Number(RealPart);
					}
				}

				this.Add = function(Complex) {
					return Mathematics.Complex(RealPart + Complex.RealPart(), ImaginaryPart + Complex.ImaginaryPart());
				};

				this.Subtract = function(Complex) {
					return Mathematics.Complex(RealPart - Complex.RealPart(), ImaginaryPart - Complex.ImaginaryPart());
				};

				this.Multiply = function(Complex) {
					return Mathematics.Complex(RealPart * Complex.RealPart() - ImaginaryPart * Complex.ImaginaryPart(), RealPart * Complex.ImaginaryPart() + ImaginaryPart * Complex.RealPart());
				};

				this.Divide = function(Complex) {
					return Mathematics.Complex((RealPart * Complex.RealPart() + ImaginaryPart * Complex.ImaginaryPart()) / (Mathematics.Square(Complex.RealPart()) + Mathematics.Square(Complex.ImaginaryPart())), (ImaginaryPart * Complex.RealPart() - RealPart * Complex.ImaginaryPart()) / (Mathematics.Square(Complex.RealPart()) + Mathematics.Square(Complex.ImaginaryPart())));
				};

				this.Conjugate = function() {
					return Mathematics.Complex(RealPart, -ImaginaryPart);
				};

				this.Stringify = function() {
					return RealPart.toString() + ' + ' + ImaginaryPart.toString() + 'i';
				}
			};

			return new Complex(A, B);
		},

		// Linear transformation
		RotateTransform: function(Angle, Point) {
			return Mathematics.Matrix([[Mathematics.Cos(Angle), -Mathematics.Sin(Angle)], [Mathematics.Sin(Angle), Mathematics.Cos(Angle)]]).Multiply(Mathematics.Matrix([[Point[0]], [Point[1]]])).ColumnVector(1).Serialize();
		},

		// Permutation and combination
		Permutation: function(N, R) {
			return Mathematics.Divide(Mathematics.Factorial(N), Mathematics.Factorial(Mathematics.Subtract(N, R)));
		},
		RepeatedPermutation: function(N, R) {
			return Mathematics.Power(N, R);
		},
		Combination: function(N, R) {
			return Mathematics.Divide(Mathematics.Permutation(N, R), Mathematics.Factorial(R));
		},
		RepeatedCombination: function(N, R) {
			return Mathematics.Combination(Mathematics.Subtract(Mathematics.Add(N, R), 1), R);
		}
	};

	// Mathematical constants
	Mathematics.E = Mathematics.Sum(0, Infinity, function(N) { return 1 / Mathematics.Factorial(N) });
	Mathematics.PI = Mathematics.Sum(0, Infinity, function(K) { return Mathematics.Power(16, -K) * (4 / (8 * K + 1) - 2 / (8 * K + 4) - 1 / (8 * K + 5) - 1 / (8 * K + 6)) });
	Mathematics.PHI = (1 + Mathematics.Root(5)) / 2;

	// (Extend) Linear algebra - Matrix
	Mathematics.Matrix.Identity = function(Size) {
		var Values = [];

		for(var RowIterator = 0; RowIterator < Size; RowIterator++) {
			Values[RowIterator] = [];

			for(var ColumnIterator = 0; ColumnIterator < Size; ColumnIterator++) {
				Values[RowIterator][ColumnIterator] = RowIterator == ColumnIterator ? 1 : 0;
			}
		}

		return Mathematics.Matrix(Values);
	};

	Mathematics.Matrix.Zero = function(Size) {
		var Values = [];

		for(var RowIterator = 0; RowIterator < Size; RowIterator++) {
			Values[RowIterator] = [];

			for(var ColumnIterator = 0; ColumnIterator < Size; ColumnIterator++) {
				Values[RowIterator][ColumnIterator] = 0;
			}
		}

		return Mathematics.Matrix(Values);
	};

	// Extend to default classes
	Number.prototype.Stringify = function() {
		return this.toString();
	}

	window.Mathematics = Mathematics;
})();