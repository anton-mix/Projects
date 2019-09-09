using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Shapes;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;

namespace prob1 {
	enum LineDir {
		LineT = 0,
		LineTR = 1,
		LineR = 2,
		LineBR = 3,
		LineB = 4,
		LineBL = 5,
		LineL = 6,
		LineTL = 7,
	}
	class Node {
		Point center;
		int rad = 0, scale = 0, countFree = 7;
		double depth = 0;
		String Name = "Привет мир";
		TextBlock NameText = new TextBlock();
		Ellipse CenterF = new Ellipse();
		public Line[] LineF = new Line[8], NameLine = new Line[2];
		Grid Perent;
		int[] freeL = new int[8];
		public int[] LengthL = new int[8];
		public Point[] endLines = new Point[8];
		TranslateTransform Trans;
		public Node(int X,int Y,int Rad,int Scale,Grid Parent,double Depth,TranslateTransform Trans) {
			center = new Point(X,Y);
			this.rad = Rad;
			this.scale = Scale;
			this.depth = Depth;
			this.Perent = Parent;
			this.Trans = Trans;
			CenterF.Fill = CenterF.Stroke = Brushes.Black;
			CenterF.StrokeThickness = 1;
			CenterF.Width = CenterF.Height = Scale * Rad * 2;
			CenterF.Margin = new Thickness(X - rad,Y - rad,0,0);
			CenterF.HorizontalAlignment = HorizontalAlignment.Left;
			CenterF.VerticalAlignment = VerticalAlignment.Top;
			CenterF.RenderTransform = Trans;
			Perent.Children.Add(CenterF);
			DrawName();
			for (int i = 0; i < 8; i++) {
				LengthL[i] = 0;
				freeL[i] = i;
			}
		}
		public int CountFree {
			get { return countFree; }
			set { countFree = value; }
		}
		public int Rad {
			get { return rad; }
			set { if (value < 0)
					rad = 0;
				else
					rad = value;
			}
		}
		public int Scale {
			get { return scale; }
			set {
				if (value < 0)
					scale = 0;
				else
					scale = value;
			}
		}
		public double Depth {
			get { return depth; }
			set { if (value < 0)
					depth = 0;
				else
					depth = value;
				for (int i = 0; i < 8; i++) {
					if (LineF[i] != null)
						LineF[i].StrokeThickness = depth;
				}
			}
		}
		public Point Center {
			get { return center; }
			set {
				CenterF.Margin = new Thickness(value.X - rad,value.Y - rad,0,0);
				for(int i = 0; i < 8; i++) {
					if(LineF[i] != null) {
						LineF[i].X1 = value.X; LineF[i].Y1 = value.Y;
						switch ((LineDir)i) {
							case LineDir.LineT:
								LineF[i].X2 = LineF[i].X2 - center.X + value.X;
								LineF[i].Y2 = LineF[i].Y2 - center.Y + value.Y;
								endLines[i] = new Point(LineF[i].X2,LineF[i].Y2);
								break;
							case LineDir.LineL:
								LineF[i].X2 = LineF[i].X2 - center.X + value.X;
								LineF[i].Y2 = LineF[i].Y2 - center.Y + value.Y;
								endLines[i] = new Point(LineF[i].X2,LineF[i].Y2);
								break;
							case LineDir.LineR:
								LineF[i].X2 = LineF[i].X2 - center.X + value.X;
								LineF[i].Y2 = LineF[i].Y2 - center.Y + value.Y;
								endLines[i] = new Point(LineF[i].X2,LineF[i].Y2);
								break;
							case LineDir.LineB:
								LineF[i].X2 = LineF[i].X2 - center.X + value.X;
								LineF[i].Y2 = LineF[i].Y2 - center.Y + value.Y;
								endLines[i] = new Point(LineF[i].X2,LineF[i].Y2);
								break;
							case LineDir.LineTL:
								LineF[i].X2 = LineF[i].X2 - center.X + value.X;
								LineF[i].Y2 = LineF[i].Y2 - center.Y + value.Y;
								endLines[i] = new Point(LineF[i].X2,LineF[i].Y2);
								break;
							case LineDir.LineTR:
								LineF[i].X2 = LineF[i].X2 - center.X + value.X;
								LineF[i].Y2 = LineF[i].Y2 - center.Y + value.Y;
								endLines[i] = new Point(LineF[i].X2,LineF[i].Y2);
								break;
							case LineDir.LineBL:
								LineF[i].X2 = LineF[i].X2 - center.X + value.X;
								LineF[i].Y2 = LineF[i].Y2 - center.Y + value.Y;
								endLines[i] = new Point(LineF[i].X2,LineF[i].Y2);
								break;
							case LineDir.LineBR:
								LineF[i].X2 = LineF[i].X2 - center.X + value.X;
								LineF[i].Y2 = LineF[i].Y2 - center.Y + value.Y;
								endLines[i] = new Point(LineF[i].X2,LineF[i].Y2);
								break;
						}
					}
				}
				center = value;
			}
		}
		public int NextLine(int rand) {
			int buff = freeL[rand];
			freeL[rand] = freeL[countFree];
			countFree--;
			return buff;
		}
		public void DrawLine(LineDir NomLine,int Length,int offsetX, int offsetY) {
			int BuffDir = (int)NomLine;
			LengthL[BuffDir] = Length;
			if (LineF[BuffDir] == null) {
				LineF[BuffDir] = new Line();
				LineF[BuffDir].X1 = center.X;
				LineF[BuffDir].Y1 = center.Y;
				LineF[BuffDir].StrokeThickness = depth;
				LineF[BuffDir].Stroke = Brushes.Black;
				LineF[BuffDir].HorizontalAlignment = HorizontalAlignment.Left;
				LineF[BuffDir].VerticalAlignment = VerticalAlignment.Top;
				LineF[BuffDir].RenderTransform = Trans;
				Perent.Children.Add(LineF[BuffDir]);
			}
			DrawName();
			switch (NomLine) {
				case LineDir.LineT:
					LineF[BuffDir].X2 = center.X + offsetX;
					LineF[BuffDir].Y2 = center.Y - Length + offsetY;
					endLines[BuffDir] = new Point(LineF[BuffDir].X2,LineF[BuffDir].Y2);
					break;
				case LineDir.LineL:
					LineF[BuffDir].X2 = center.X - Length + offsetX;
					LineF[BuffDir].Y2 = center.Y + offsetY;
					endLines[BuffDir] = new Point(LineF[BuffDir].X2,LineF[BuffDir].Y2);
					break;
				case LineDir.LineR:
					LineF[BuffDir].X2 = center.X + Length + offsetX;
					LineF[BuffDir].Y2 = center.Y + offsetY;
					endLines[BuffDir] = new Point(LineF[BuffDir].X2,LineF[BuffDir].Y2);
					break;
				case LineDir.LineB:
					LineF[BuffDir].X2 = center.X + offsetX;
					LineF[BuffDir].Y2 = center.Y + Length + offsetY;
					endLines[BuffDir] = new Point(LineF[BuffDir].X2,LineF[BuffDir].Y2);
					break;
				case LineDir.LineTL:
					LineF[BuffDir].X2 = center.X - Length + offsetX;
					LineF[BuffDir].Y2 = center.Y - Length + offsetY;
					endLines[BuffDir] = new Point(LineF[BuffDir].X2,LineF[BuffDir].Y2);
					break;
				case LineDir.LineTR:
					LineF[BuffDir].X2 = center.X + Length + offsetX;
					LineF[BuffDir].Y2 = center.Y - Length + offsetY;
					endLines[BuffDir] = new Point(LineF[BuffDir].X2,LineF[BuffDir].Y2);
					break;
				case LineDir.LineBL:
					LineF[BuffDir].X2 = center.X - Length + offsetX;
					LineF[BuffDir].Y2 = center.Y + Length + offsetY;
					endLines[BuffDir] = new Point(LineF[BuffDir].X2,LineF[BuffDir].Y2);
					break;
				case LineDir.LineBR:
					LineF[BuffDir].X2 = center.X + Length + offsetX;
					LineF[BuffDir].Y2 = center.Y + Length + offsetY;
					endLines[BuffDir] = new Point(LineF[BuffDir].X2,LineF[BuffDir].Y2);
					break;
			}
		}
		public void MoveNode(double offsetX,double offsetY) {
			center.X = center.X + offsetX;
			center.Y = center.Y + offsetY;
			CenterF.Margin = new Thickness(CenterF.Margin.Left + offsetX,CenterF.Margin.Top+offsetY,0,0);
			DrawName();
			for (int i = 0; i < 8; i++) {
				if(LineF[i] != null) {
					if (LengthL[i] == -2) {
						LineF[i].X2 = center.X;
						LineF[i].Y2 = center.Y;
					}
					else {
						LineF[i].X1 = center.X;
						LineF[i].Y1 = center.Y;
					}
				}
			}
		}
		void DrawName() {
			if (NameLine[0] == null) {
				for (int i = 0; i < 2; i++) {
					NameLine[i] = new Line();
					NameLine[i].StrokeThickness = depth;
					NameLine[i].Stroke = Brushes.Black;
					NameLine[i].HorizontalAlignment = HorizontalAlignment.Left;
					NameLine[i].VerticalAlignment = VerticalAlignment.Top;
					NameLine[i].RenderTransform = Trans;
					Perent.Children.Add(NameLine[i]);
				}
				NameText.Text = Name;
				NameText.HorizontalAlignment = HorizontalAlignment.Left;
				NameText.VerticalAlignment = VerticalAlignment.Top;
				NameText.RenderTransform = Trans;
				NameText.Measure(new Size(Double.PositiveInfinity,Double.PositiveInfinity));
				NameText.Arrange(new Rect(NameText.DesiredSize));
				Perent.Children.Add(NameText);
			}
			/*LineT = 0,
		LineTR = 1,
		LineR = 2,
		LineBR = 3,
		LineB = 4,
		LineBL = 5,
		LineL = 6,
		LineTL = 7,*/
			NameLine[0].X1 = center.X;
			NameLine[0].Y1 = center.Y;
			if (LineF[1] == null) {
				NameText.Margin = new Thickness(center.X + 35,center.Y - 20 - NameText.ActualHeight,0,0);
				NameLine[0].X2 = NameLine[1].X1 = center.X + 20;
				NameLine[0].Y2 = NameLine[1].Y1 = NameLine[1].Y2 = center.Y - 20;
				NameLine[1].X2 = center.X + 35 + NameText.ActualWidth;
			}
			else if (LineF[7] == null) {
				NameText.Margin = new Thickness(center.X - 35 - NameText.ActualWidth,center.Y - 20 - NameText.ActualHeight,0,0);
				NameLine[0].X2 = NameLine[1].X1 = center.X - 20;
				NameLine[0].Y2 = NameLine[1].Y1 = NameLine[1].Y2 = center.Y - 20;
				NameLine[1].X2 = center.X - 35 - NameText.ActualWidth;
			}
			else if (LineF[3] == null) {
				NameText.Margin = new Thickness(center.X + 35,center.Y + 20 - NameText.ActualHeight,0,0);
				NameLine[0].X2 = NameLine[1].X1 = center.X + 20;
				NameLine[0].Y2 = NameLine[1].Y1 = NameLine[1].Y2 = center.Y + 20;
				NameLine[1].X2 = center.X + 35 + NameText.ActualWidth;
			}
			else if(LineF[5] == null) {
				NameText.Margin = new Thickness(center.X - 35 - NameText.ActualWidth,center.Y + 20 - NameText.ActualHeight,0,0);
				NameLine[0].X2 = NameLine[1].X1 = center.X - 20;
				NameLine[0].Y2 = NameLine[1].Y1 = NameLine[1].Y2 = center.Y + 20;
				NameLine[1].X2 = center.X - 35 - NameText.ActualWidth;
			}
			else {
				NameText.Margin = new Thickness(center.X + 35,center.Y - 20 - NameText.ActualHeight,0,0);
				NameLine[0].X1 = center.X;
				NameLine[0].Y1 = center.Y;
				NameLine[0].X2 = NameLine[1].X1 = center.X + 20;
				NameLine[0].Y2 = NameLine[1].Y1 = center.Y - 20;
				NameLine[1].X2 = center.X + 35 + NameText.ActualWidth;
				NameLine[1].Y2 = center.Y - 20;
			}
		}
	}
}
