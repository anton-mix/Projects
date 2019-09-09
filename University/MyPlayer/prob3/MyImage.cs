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
	class MyImage {
		double width, height,pictureWidth, picturHeight;
		string nameImage;
		Rect border;
		DrawingGroup image;
		BitmapImage picture;
		FormattedText TextF;
		int Point = 14;
		public MyImage(BitmapImage Bimage,Rect border, DrawingGroup MainDraw) {
			this.picture = Bimage;
			this.border = border;
			image = new DrawingGroup();

			MainDraw.Children.Add(image);

			Update();
		}

		void Update() {
			image.Children.Clear();
			if (picture == null) {
				TextF = new FormattedText("Not Image", CultureInfo.GetCultureInfo("en-us"), FlowDirection.LeftToRight, new Typeface("Areal"), Point, Brushes.Black);
				Geometry TextG = TextF.BuildGeometry(new Point(border.X + (border.Width / 2 - TextF.Width/2),border.Y + (border.Height / 2 - TextF.Height / 2)));
				GeometryDrawing TextGD = new GeometryDrawing(Brushes.Black, new Pen(Brushes.Black, 0), TextG);
				GeometryDrawing Rectangle = new GeometryDrawing(Brushes.Transparent, new Pen(Brushes.Black, 1), new RectangleGeometry(border));
				image.Children.Add(Rectangle);
				image.Children.Add(TextGD);
				pictureWidth = border.Width;
				picturHeight = border.Height;
			}
			else {
				image.Children.Add(new ImageDrawing(picture, border));
				pictureWidth = picture.Width;
				picturHeight = picture.Height;
			}
			width = border.Width;
			height = border.Height;
		}

		public BitmapImage Picture {
			get { return picture; }
			set {
				picture = value;
				Update();
			}
		}

		public double Width {
			get { return width; }
		}

		public double Height {
			get { return height; }
		}

		public DrawingGroup Image {
			get { return image; }
		}

		public Rect Border {
			get { return border; }
			set {
				border = value;
				Update();
			}
		}

	}
}
