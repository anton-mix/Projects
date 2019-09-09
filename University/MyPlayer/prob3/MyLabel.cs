using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Media.Imaging;
using System.Diagnostics;
using System.Windows.Automation;
using System.Windows.Controls;
using System.Windows.Controls.Primitives;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Ink;
using System.Windows.Input;
using System.Windows.Markup;
using System.Windows.Media;
using System.Windows.Media.Animation;
using System.Windows.Media.Effects;
using System.Windows.Media.Media3D;
using System.Windows.Media.TextFormatting;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Windows.Shell;
using System.Globalization;

namespace prob3 {
	class MyLabel {
		string text, header;
		float div = 0.3F;
		double width, height, div2 = 5, div3 = 5,widthBody, widthHeader;
		FormattedText TextF, HeaderF;
		DrawingGroup label;
		int Point;
		Point coord;
		Color colorHeadr;
		public MyLabel(string header, string text, int Point,double X, double Y, DrawingGroup MainDraw) {
			this.text = text;
			this.header = header;
			this.Point = Point;
			this.label = new DrawingGroup();
			MainDraw.Children.Add(label);
			coord = new Point(X,Y);
			colorHeadr = Color.FromArgb(255, 0, 147, 221);
			widthBody = -1;
			widthHeader = -1;
			Update(); 
		}
		void Update() {
			label.Children.Clear();
			TextF = new FormattedText(text, CultureInfo.GetCultureInfo("en-us"), FlowDirection.LeftToRight, new Typeface("Areal"), Point, Brushes.Black);
			HeaderF = new FormattedText(header, CultureInfo.GetCultureInfo("en-us"), FlowDirection.LeftToRight, new Typeface("Areal"), Point, Brushes.Black);

			if (widthBody != -1) {
				while (true) {
					if (TextF.Width > widthBody - 4) {
						text = text.Remove(text.Length - 4) + "...";
						TextF = new FormattedText(text, CultureInfo.GetCultureInfo("en-us"), FlowDirection.LeftToRight, new Typeface("Areal"), Point, Brushes.Black);
					}
					else {
						break;
					}
				}
			}

			Geometry HeaderG = HeaderF.BuildGeometry(new Point(coord.X + div3, coord.Y + div3));
			Geometry TextG;
			if (widthHeader == -1) {
				TextG = TextF.BuildGeometry(new Point(coord.X + div3 * 3 + div2 + HeaderF.Width, coord.Y + div3));
			}
			else {
				TextG = TextF.BuildGeometry(new Point(coord.X + div3 * 1 + div2 + widthHeader, coord.Y + div3));
			}
			GeometryDrawing TextGD = new GeometryDrawing(Brushes.Black, new Pen(Brushes.Black, 0), TextG);
			GeometryDrawing HeaderGD = new GeometryDrawing(Brushes.Black, new Pen(Brushes.Black, 0), HeaderG);

			PathSegmentCollection myPathSegmentCollection = new PathSegmentCollection();
			if (widthHeader == -1) {
				myPathSegmentCollection.Add(new LineSegment(new Point(coord.X, coord.Y + (HeaderF.Height + div3 * 2) * (1 - div)), false));
				myPathSegmentCollection.Add(new LineSegment(new Point(coord.X + (HeaderF.Height + div3 * 2) * div, coord.Y + (HeaderF.Height + div3 * 2)), false));
				myPathSegmentCollection.Add(new LineSegment(new Point(coord.X + HeaderF.Width + div3 * 2, coord.Y + HeaderF.Height + div3 * 2), false));
				myPathSegmentCollection.Add(new LineSegment(new Point(coord.X + HeaderF.Width + div3 * 2, coord.Y), false));
			}
			else {
				myPathSegmentCollection.Add(new LineSegment(new Point(coord.X, coord.Y + (HeaderF.Height + div3 * 2) * (1 - div)), false));
				myPathSegmentCollection.Add(new LineSegment(new Point(coord.X + (HeaderF.Height + div3 * 2) * div, coord.Y + (HeaderF.Height + div3 * 2)), false));
				myPathSegmentCollection.Add(new LineSegment(new Point(coord.X + widthHeader, coord.Y + HeaderF.Height + div3 * 2), false));
				myPathSegmentCollection.Add(new LineSegment(new Point(coord.X + widthHeader, coord.Y), false));
			}

			PathFigureCollection myPathFigureCollection = new PathFigureCollection();
			myPathFigureCollection.Add(new PathFigure(new Point(coord.X, coord.Y), myPathSegmentCollection, false));

			GeometryDrawing Path1 = new GeometryDrawing(new SolidColorBrush(colorHeadr), new Pen(Brushes.Black, 1), new PathGeometry(myPathFigureCollection));
			GeometryDrawing Rectangle1;
			if (widthBody == -1) {
				if (widthHeader == -1) {
					Rectangle1 = new GeometryDrawing(Brushes.Transparent, new Pen(Brushes.Black, 1), new RectangleGeometry(new Rect(coord.X + div3 * 2 + div2 + HeaderF.Width, coord.Y, TextF.Width + div3 * 2, TextF.Height + div3 * 2)));
				}
				else {
					Rectangle1 = new GeometryDrawing(Brushes.Transparent, new Pen(Brushes.Black, 1), new RectangleGeometry(new Rect(coord.X + div2 + widthHeader, coord.Y, TextF.Width + div3 * 2, TextF.Height + div3 * 2)));
				}
			}
			else {
				if (widthHeader == -1) {
					Rectangle1 = new GeometryDrawing(Brushes.Transparent, new Pen(Brushes.Black, 1), new RectangleGeometry(new Rect(coord.X + div3 * 2 + div2 + HeaderF.Width, coord.Y, widthBody, TextF.Height + div3 * 2)));
				}
				else {
					Rectangle1 = new GeometryDrawing(Brushes.Transparent, new Pen(Brushes.Black, 1), new RectangleGeometry(new Rect(coord.X + div2 + widthHeader, coord.Y, widthBody, TextF.Height + div3 * 2)));
				}
			}

			label.Children.Add(Rectangle1);
			label.Children.Add(Path1);
			label.Children.Add(TextGD);
			label.Children.Add(HeaderGD);
			width = HeaderF.Width + TextF.Width + div2 + div3 * 4;
			height = ((HeaderF.Height > TextF.Height) ? HeaderF.Height : TextF.Height) + div3 * 2;
		}
		public string Text {
			get { return text; }
			set {
				text = value;
				Update();
			}
		}

		public string Header {
			get { return header; }
			set {
				header = value;
				Update();
			}
		}
		public DrawingGroup Label {
			get { return label; }
		}

		public double Width {
			get { return width; }
		}

		public double Height {
			get { return height; }
		}

		public Color ColorHeadr {
			get { return colorHeadr; }
			set {
				colorHeadr = value;
				Update();
			}
		}

		public double WidthHeader {
			get { return widthHeader; }
			set {
				widthHeader = value;
				Update();
			}
		}

		public double WidthBody {
			get { return widthBody; }
			set {
				widthBody = value;
				Update();
			}
		}
	}
}
