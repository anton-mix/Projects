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
	class MyButton {
		BitmapImage buttonOn, buttonOff;
		DrawingGroup button;
		string text;
		Rect border;
		double height, width;
		int point = 15;
		bool stateButton;
		public delegate void Event();
		public Event EventButtonKeyDown;
		public Event EventButtonKeyUp;
		public Event EventButtonKeyMove;
		bool enabled;

		public MyButton(Rect Border,string text,int Point, DrawingGroup MainDraw) {
			this.border = Border;
			this.text = text;
			this.point = Point;
			button = new DrawingGroup();
			MainDraw.Children.Add(button);
			buttonOff = new BitmapImage();
			buttonOff.BeginInit();
			buttonOff.UriSource = new Uri("image/button1_2.png", UriKind.Relative);
			buttonOff.EndInit();
			
			buttonOn = new BitmapImage();
			buttonOn.BeginInit();
			buttonOn.UriSource = new Uri("image/button2_2.png", UriKind.Relative);
			buttonOn.EndInit();

			button.Opacity = 1;
			stateButton = false;
			enabled = true;
			Update(buttonOff);
		}

		public MyButton(Rect Border, string text, int Point, BitmapImage buttonOn, BitmapImage buttonOff) {
			this.border = Border;
			this.text = text;
			this.point = Point;
			button = new DrawingGroup();

			this.buttonOff = buttonOff;
			this.buttonOn = buttonOn;

			Update(buttonOff);
		}

		public void Button_MouseKeyDown(Point point) {
			if (enabled && Test(border,point)) {
				stateButton = true;
				Update(buttonOn);
				if(EventButtonKeyDown != null)
					EventButtonKeyDown();
			}
		}

		public void Button_MouseKeyUp(Point point) {
			if (enabled && stateButton && Test(border, point)) {
				stateButton = false;
				Update(buttonOff);
				if (EventButtonKeyUp != null)
					EventButtonKeyUp();
			}
		}

		public void Button_MouseKeyMove(Point point) {
			if (enabled && stateButton && !Test(border, point)) {
				stateButton = false;
				Update(buttonOff);
				if (EventButtonKeyMove != null)
					EventButtonKeyMove();
			}
		} 

		bool Test(Rect Border, Point Key) {
			if (Key.X < Border.BottomRight.X && Key.X > Border.TopLeft.X && Key.Y < Border.BottomRight.Y && Key.Y > Border.TopLeft.Y) {
				return true;
			}
			return false;
		}

		void Update(BitmapImage ImageButton) {
			FormattedText TextF = new FormattedText(text, CultureInfo.GetCultureInfo("en-us"), FlowDirection.LeftToRight, new Typeface("Areal"), point, Brushes.Black);

			Geometry TextG = TextF.BuildGeometry(new Point(border.X + (border.Width / 2 - TextF.Width / 2), border.Y + (border.Height / 2 - TextF.Height / 2)));

			GeometryDrawing TextGD = new GeometryDrawing(Brushes.White, new Pen(Brushes.White, 0), TextG);

			button.Children.Clear();
			button.Children.Add(new ImageDrawing(ImageButton, border));
			button.Children.Add(TextGD);
			width = border.Width;
			height = border.Height;
		}

		public double Width {
			get { return width; }
		}

		public double Height {
			get { return height; }
		}

		public DrawingGroup Button {
			get { return button; }
		}

		public Rect Border {
			get { return border; }
			set {
				border = value;
				Update(buttonOff);
			}
		}

		public int Point {
			get { return point; }
			set {
				point = value;
				Update(buttonOff);
			}
		}
		public string Text {
			get { return text; }
			set {
				text = value;
				Update(buttonOff);
			}
		}

		public bool Enabled {
			get { return enabled; }
			set {
				enabled = value;
				button.Opacity = (enabled) ? 1 : 0.8;
			}
		}

	}
}
