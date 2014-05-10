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
		Power: Math.pow,
		Root: function(Value, Powered) {
			if(!Powered) {
				Powered = 2;
			}

			return Mathematics.Power(Value, 1 / Powered);
		},
		Square: function(X) {
			return Mathematics.Power(X, 2);
		},
		Cube: function(X) {
			return Mathematics.Power(X, 3);
		},
		Sin: Math.sin,
		Cos: Math.cos,
		Tan: Math.tan,
		Sec: function(X) {
			return 1 / Math.cos(X);
		},
		Csc: function(X) {
			return 1 / Math.sin(X);	
		},
		Cot: function(X) {
			return 1 / Math.tan(X);
		},
		Arcsin: Math.asin,
		Arccos: Math.acos,
		Arctan: Math.atan,
		Arcsec: function(X) {
			return 1 / Math.acos(X);
		},
		Arccsc: function(X) {
			return 1 / Math.asin(X);	
		},
		Arccot: function(X) {
			return 1 / Math.atan(X);
		},
		Sinh: function(X) {
			return (Mathematics.Exp(X) - Mathematics.Exp(-X)) / 2;
		},
		Cosh: function(X) {
			return (Mathematics.Exp(X) + Mathematics.Exp(-X)) / 2;
		},
		Tanh: function(X) {
			return Mathematics.Sinh(X) / Mathematics.Cosh(X);
		},
		Sech: function(X) {
			return 1 / Mathematics.Cosh(X);
		},
		Csch: function(X) {
			return 1 / Mathematics.Sinh(X);	
		},
		Coth: function(X) {
			return 1 / Mathematics.Tanh(X);
		},
		Arcsinh: function(X) {
			return Mathematics.Ln(X + Mathematics.Root(Mathematics.Square(X) + 1));
		},
		Arccosh: function(X) {
			return X >= 1 ? Mathematics.Ln(X + Mathematics.Root(Mathematics.Square(X) - 1)) : undefined;
		},
		Arctanh: function(X) {
			return Mathematics.Abs(X) < 1 ? Mathematics.Ln((X + 1) / (1 - X)) / 2 : undefined;
		},
		Arcsech: function(X) {
			return 0 < X && X <= 1 ? Mathematics.Ln(1 / X + Mathematics.Root(1 - Mathematics.Square(X)) / X) : undefined;
		},
		Arccsch: function(X) {
			return X != 0 ? Mathematics.Ln(1 / X + Mathematics.Root(1 + Mathematics.Square(X)) / Mathematics.Abs(X)) : undefined;
		},
		Arccoth: function(X) {
			if(Mathematics.Abs(X) > 1) {
				return Mathematics.Ln((X + 1) / (X - 1)) / 2;
			}
		},
		Exp: Math.exp,
		Log: function(X) {
			return Math.log(X) / Math.log(10);
		},
		Ln: Math.log,
		Abs: Math.abs,
		Sign: function(X) {
			return X == 0 ? X : X / Mathematics.Abs(X);
		},
		Floor: Math.floor,
		Ceil: Math.ceil,
		Sum: function(From, To, Fx) {
			var Summation = 0;

			for(var Iterator = From; Iterator <= To; Iterator++) {
				if(Mathematics.Abs(To) == Infinity) {
					if(Summation != Summation + Fx(Iterator)) {
						Summation += Fx(Iterator);
					} else {
						break;
					}
				} else {
					Summation += Fx(Iterator);	
				}
			}

			return Summation;
		},
		Product: function(From, To, Fx) {
			var Production = 1;

			for(var Iterator = From; Iterator <= To; Iterator++) {
				if(Mathematics.Abs(To) == Infinity) {
					if(Production != Production * Fx(Iterator)) {
						Production *= Fx(Iterator);
					} else {
						break;
					}
				} else {
					Production *= Fx(Iterator);
				}
			}

			return Production;
		},
		Zeta: function(S) {
			if(S > 1) {
				return S == Infinity ? 1 : Mathematics.Sum(1, Infinity, function(N) { return 1 / Mathematics.Power(N, S) });	
			} else {
				return undefined;
			}
		},
		Gamma: function(T) {
			if(T > 0 || T < 0 && !Mathematics.isInteger(T)) {
				return Mathematics.isInteger(T) ? Mathematics.Product(1, T - 1, function(N) { return N }) : Mathematics.Product(1, Infinity, function(N) { return Mathematics.Power(1 + 1 / N, T) / (1 + T / N) }) / T;	
			} else {
				return undefined;
			}
		},
		Factorial: function(N) {
			return Mathematics.isInteger(N) ? Mathematics.Gamma(N + 1) : undefined;
		},
		isInteger: function(X) {
			return X % 1 == 0;
		},
		Radian: function(Degree) {
			return Degree / 180 * Mathematics.PI;
		},
		Degree: function(Radian) {
			return Radian * 180 / Mathematics.PI;
		},
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
				}

				this.ColumnVector = function(Column) {
					var Vector = [];

					for(var RowIterator = 1; RowIterator <= Rows; RowIterator++) {
						Vector.push([this.Item(RowIterator, Column)]);
					}

					return Mathematics.Matrix(Vector);
				}

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

					return Mathematics.Matrix(Minor);
				}

				this.Cofactor = function(Row, Column) {
					return Mathematics.Power(-1, Row + Column) * this.Minor(Row, Column).Determinant();
				}

				this.Adjugate = function() {
					var Matrix = [];

					for(var RowIterator = 1; RowIterator <= Rows; RowIterator++) {
						Matrix.push([]);

						for(var ColumnIterator = 1; ColumnIterator <= Columns; ColumnIterator++) {
							Matrix[Matrix.length - 1].push(this.Cofactor(RowIterator, ColumnIterator));
						}
					}

					return Mathematics.Matrix(Matrix).Transpose();
				}

				this.ScalarMultiply = function(Scalar) {
					for(var RowIterator = 0; RowIterator < Rows; RowIterator++) {
						for(var ColumnIterator = 0; ColumnIterator < Columns; ColumnIterator++) {
							ValueSet[RowIterator][ColumnIterator] *= Scalar;
						}
					}
				}

				this.Transpose = function() {
					var Values = [];

					for(var RowIterator = 0; RowIterator < Columns; RowIterator++) {
						Values[RowIterator] = [];

						for(var ColumnIterator = 0; ColumnIterator < Rows; ColumnIterator++) {
							Values[RowIterator][ColumnIterator] = ValueSet[ColumnIterator][RowIterator];
						}
					}

					return Mathematics.Matrix(Values);
				}

				this.Multiply = function(Matrix) {
					if(Columns == Matrix.Rows) {
						var Induction = [], Transposal = Matrix.Transpose();

						for(var RowIterator = 0; RowIterator < Rows; RowIterator++) {
							Induction[RowIterator] = [];

							for(var ColumnIterator = 0; ColumnIterator < Matrix.Columns; ColumnIterator++) {
								Induction[RowIterator][ColumnIterator] = 0;

								for(var Iterator = 1; Iterator <= Columns; Iterator++) {
									Induction[RowIterator][ColumnIterator] += this.Item(RowIterator + 1, Iterator) * Transposal.Item(ColumnIterator + 1, Iterator);
								}
							}
						}

						return Mathematics.Matrix(Induction);
					}
				}

				if(Columns == Rows) {
					this.Trace = function(Matrix) {
						var Trace = 0;

						for(var Iterator = 1; Iterator <= Rows; Iterator++) {
							Trace += this.Item(Iterator, Iterator);
						}

						return Trace;
					};

					this.Determinant = function(Matrix) {
						if(Rows == 1) {
							return this.Item(1, 1);
						} else {
							var Production = 1, Difference = 1;

							for(var ColumnIterator = 1; ColumnIterator <= Columns; ColumnIterator++) {
								Production *= this.Item(ColumnIterator, ColumnIterator);
								Difference *= this.Item(ColumnIterator, Columns + 1 - ColumnIterator);
							}

							return Production - Difference;
						}
					};

					if(this.Determinant() != 0) {
						this.InverseMatrix = function() {
							var InverseMatrix = this.Adjugate();
							InverseMatrix.ScalarMultiply(1 / this.Determinant());

							return InverseMatrix;
						}
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
					}
				}

				this.Stringify = function() {
					var MatrixString = '[';

					for(var RowIterator = 1; RowIterator <= Rows; RowIterator++) {
						MatrixString += RowIterator == 1 ? '' : ' ';

						for(var ColumnIterator = 1; ColumnIterator <= Columns; ColumnIterator++) {
							MatrixString += this.Item(RowIterator, ColumnIterator).toString() + (RowIterator == Rows && ColumnIterator == Columns ? '' : ', ');
						}

						MatrixString += RowIterator == Rows ? ']' : '\n';
					}

					console.log(MatrixString);
				}
			};

			return new Matrix(Values);
		},
		RotateTransform: function(Angle, Point) {
			return Mathematics.Matrix([[Mathematics.Cos(Angle), -Mathematics.Sin(Angle)], [Mathematics.Sin(Angle), Mathematics.Cos(Angle)]]).Multiply(Mathematics.Matrix([[Point[0]], [Point[1]]])).ColumnVector(1).Serialize();
		},
		Permutate: function(N, R) {
			return Mathematics.Factorial(N) / Mathematics.Factorial(N - R);
		},
		Combinate: function(N, R) {
			return Mathematics.Permutate(N, R) / Mathematics.Factorial(R);
		}
	};

	Mathematics.E = Mathematics.Sum(0, Infinity, function(N) { return 1 / Mathematics.Factorial(N) });
	Mathematics.PI = Mathematics.Sum(0, Infinity, function(K) { return Mathematics.Power(16, -K) * (4 / (8 * K + 1) - 2 / (8 * K + 4) - 1 / (8 * K + 5) - 1 / (8 * K + 6)) });
	Mathematics.PHI = (1 + Mathematics.Root(5)) / 2;

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

	window.Mathematics = Mathematics;
})();