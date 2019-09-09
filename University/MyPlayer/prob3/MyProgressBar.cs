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
	class MyProgressBar {
		enum RectBorder : int {
			Form = 0,
			Ellement = 1
		}


		BitmapImage BarImage, ItemImage;
		DrawingGroup Bar;
		string text;
		double height, width, progress, min, max, barOpacity;
		Point coord;
		Rect[] borders;
		public MyProgressBar(Rect[] border, DrawingGroup MainDraw) {
			Bar = new DrawingGroup();
			MainDraw.Children.Add(Bar);
			this.borders = border;

			BarImage = new BitmapImage();
			BarImage.BeginInit();
			BarImage.UriSource = new Uri("image/ProcentBar1_2.png", UriKind.Relative);
			BarImage.EndInit();

			ItemImage = new BitmapImage();
			ItemImage.BeginInit();
			ItemImage.UriSource = new Uri("image/element3_2.png", UriKind.Relative);
			ItemImage.EndInit();

			barOpacity = 1;
			max = 100;
			min = 0;
			text = "00:01 / 02:17";

			Update();
		}

		void Update() {
			Bar.Children.Clear();
			DrawingGroup buff2 = new DrawingGroup();
			for (int i = 0; i < (int)((int)progress * 15 / (int)max); i++) {
				Rect Buff = borders[(int)RectBorder.Ellement];
				Buff.X = borders[(int)RectBorder.Form].X + 13 + i * Buff.Width * 1.1;
				Buff.Y = borders[(int)RectBorder.Form].Y + 4;
				buff2.Children.Add(new ImageDrawing(ItemImage, Buff));
				buff2.Opacity = barOpacity;
				//new Rect(100 + 7.5 + i * 22, 302, 20, 13);
			}
			Bar.Children.Add(buff2);
			Bar.Children.Add(new ImageDrawing(BarImage, borders[(int)RectBorder.Form]));

			FormattedText TextF = new FormattedText(text, CultureInfo.GetCultureInfo("en-us"), FlowDirection.LeftToRight, new Typeface("Areal"), 14, Brushes.Black);
			Geometry TextG = TextF.BuildGeometry(new Point(borders[(int)RectBorder.Form].X + (borders[(int)RectBorder.Form].Width/2 - TextF.Width / 2), borders[(int)RectBorder.Form].Y + (borders[(int)RectBorder.Form].Height / 2 - TextF.Height / 2)));
			GeometryDrawing TextGD = new GeometryDrawing(Brushes.Black, new Pen(Brushes.Black, 0), TextG);

			Bar.Children.Add(TextGD);

			height = borders[(int)RectBorder.Form].Height;
			width = borders[(int)RectBorder.Form].Width;
		}

		public double BarOpacity {
			get { return barOpacity; }
			set { barOpacity = value; }
		}

		public double Width {
			get { return width; }
		}

		public double Height {
			get { return height; }
		}

		public double Progress {
			get { return progress; }
			set {
				progress = value;
				Update();
			}
		}

		public double Max {
			get { return max; }
			set {
				max = value;
				Update();
			}
		}

		public double Min {
			get { return min; }
			set { min = value;
				Update();
			}
		}

		public string Text {
			get { return text; }
			set {
				text = value;
			}
		}
	}
}

