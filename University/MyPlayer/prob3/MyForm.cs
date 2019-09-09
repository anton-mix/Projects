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
	class MyForm {
		List<BitmapImage> Images;
		Rect[] borders;
		List<string> Paths;
		DrawingGroup form;
		string name;
		double height, width;
		Point coord;
		Image ImageForm;
		MainWindow mainWindow;
		bool stateExit, stateDrop, stateCollapse;
		public delegate void Event();
		public Event EventFormClose;
		enum Object {
			Header = 0,
			Body = 1,
			ExitI = 2,
			Collapse = 3,
			Drop = 4,
		}

		enum ArrayObjects : int {
			Header = 0,
			Body = 1,
			ExitOff = 2,
			ExitOn = 3,
			CollapseOff = 4,
			CollapseOn = 5,
			DropOff = 6,
			DropOn = 7,
		}

		public MyForm(string Name, Rect[] borders,Point coord , double width, double height , Image ImageForm, MainWindow mainWindow, DrawingGroup MainDraw) {
			this.name = Name;
			this.ImageForm = ImageForm;
			this.mainWindow = mainWindow;
			this.width = width;
			this.height = height;
			this.borders = borders;
			this.coord = coord;
			Images = new List<BitmapImage>();
			form = new DrawingGroup();
			Paths = new List<string>();
			Paths.Add("image/statusLine1.2.png");
			Paths.Add("image/background1.2.png");
			Paths.Add("image/Exit1_2.png");
			Paths.Add("image/Exit2_2.png");
			Paths.Add("image/collapse1_2.png");
			Paths.Add("image/collapse2_2.png");
			Paths.Add("image/drop1_2.png");
			Paths.Add("image/drop2_2.png");

			MainDraw.Children.Add(form);

			for (int i = 0; i < Paths.Count; i++) {
				BitmapImage BuffImage = new BitmapImage();
				BuffImage.BeginInit();
				BuffImage.UriSource = new Uri(Paths[i], UriKind.Relative);
				BuffImage.EndInit();
				Images.Add(BuffImage);
			}

			borders[(int)Object.Body].Width = Images[(int)ArrayObjects.Body].Width;
			borders[(int)Object.Body].Height = Images[(int)ArrayObjects.Body].Height;
			borders[(int)Object.ExitI].X = width - 20 - 30;
			borders[(int)Object.Collapse].X = width - 60 - 30;
			borders[(int)Object.Drop].X = width - 100 - 30;
			stateExit = stateDrop = stateCollapse = false;

			Update();
		}

		public MyForm(string Name, double width, double height, Image MainForm) {}

		public void Form_MouseDownKey(Point point) {
		
			if (Test(form.Children[(int)Object.ExitI].Bounds, point)) {
				form.Children[(int)Object.ExitI] = new ImageDrawing(Images[(int)ArrayObjects.ExitOn], borders[(int)Object.ExitI]);
				stateExit = true;
			}
			else if (Test(form.Children[(int)Object.Collapse].Bounds, point)) {
				form.Children[(int)Object.Collapse] = new ImageDrawing(Images[(int)ArrayObjects.CollapseOn], borders[(int)Object.Collapse]);
				stateCollapse = true;
			}
			else if (Test(form.Children[(int)Object.Drop].Bounds, point)) {
				form.Children[(int)Object.Drop] = new ImageDrawing(Images[(int)ArrayObjects.DropOn], borders[(int)Object.Drop]);
				stateDrop = true;
			}
			else {
				if (point.Y < 30) {
					mainWindow.DragMove();
				}
			}
		}

		public void Form__MouseUpKey(Point point) {
			if (stateExit && Test(form.Children[(int)Object.ExitI].Bounds, point)) {
				form.Children[(int)Object.ExitI] = new ImageDrawing(Images[(int)ArrayObjects.ExitOff], borders[(int)Object.ExitI]);
				stateExit = false;
				if(EventFormClose != null) {
					EventFormClose();
				}
				mainWindow.Close();
			}
			else if (stateCollapse && Test(form.Children[(int)Object.Collapse].Bounds, point)) {
				form.Children[(int)Object.Collapse] = new ImageDrawing(Images[(int)ArrayObjects.CollapseOff], borders[(int)Object.Collapse]);
				stateCollapse = false;
			}
			else if (stateDrop && Test(form.Children[(int)Object.Drop].Bounds, point)) {
				form.Children[(int)Object.Drop] = new ImageDrawing(Images[(int)ArrayObjects.DropOff], borders[(int)Object.Drop]);
				stateDrop = false;
				mainWindow.WindowState = WindowState.Minimized;
			}
		}

		public void Form__MouseMoveKey(Point point) {
			if (stateExit && !Test(form.Children[(int)Object.ExitI].Bounds, point)) {
				form.Children[(int)Object.ExitI] = new ImageDrawing(Images[(int)ArrayObjects.ExitOff], borders[(int)Object.ExitI]);
				stateExit = false;
			}
			if (stateCollapse && !Test(form.Children[(int)Object.Collapse].Bounds, point)) {
				form.Children[(int)Object.Collapse] = new ImageDrawing(Images[(int)ArrayObjects.CollapseOff], borders[(int)Object.Collapse]);
				stateCollapse = false;
			}
			if (stateDrop && !Test(form.Children[(int)Object.Drop].Bounds, point)) {
				form.Children[(int)Object.Drop] = new ImageDrawing(Images[(int)ArrayObjects.DropOff], borders[(int)Object.Drop]);
				stateDrop = false;
			}
		}

		void Update() {
			form.Children.Clear();

			double buffX = borders[(int)Object.Body].X;
			double buffY = borders[(int)Object.Body].Y;

			PathSegmentCollection myPathSegmentCollection = new PathSegmentCollection();
			myPathSegmentCollection.Add(new LineSegment(new Point(buffX + 5 + 30, buffY + 5), false));
			myPathSegmentCollection.Add(new LineSegment(new Point(buffX + 5 + 5, buffY + 5 + 30), false));
			myPathSegmentCollection.Add(new LineSegment(new Point(buffX + 5, buffY + 5 + 30), false));

			PathFigureCollection myPathFigureCollection = new PathFigureCollection();
			myPathFigureCollection.Add(new PathFigure(new Point(buffX + 5 + 25, buffY + 5), myPathSegmentCollection, false));

			GeometryDrawing Path = new GeometryDrawing(Brushes.Black, new Pen(Brushes.Black, 1), new PathGeometry(myPathFigureCollection));

			GeometryDrawing Rectangle = new GeometryDrawing(Brushes.Black, new Pen(Brushes.Black, 1), new RectangleGeometry(new Rect(buffX+ 5 + 30, buffY+5+28,borders[(int)Object.Body].Width - 60,2)));

			FormattedText TextF = new FormattedText("Media Player", CultureInfo.GetCultureInfo("en-us"), FlowDirection.LeftToRight, new Typeface("Areal"), 14, Brushes.Black);
			Geometry TextG = TextF.BuildGeometry(new Point(buffX + 5 + 30, buffY + 5 + 10));
			GeometryDrawing TextGD = new GeometryDrawing(Brushes.Black, new Pen(Brushes.Black, 0), TextG);

			form.Children.Add(new ImageDrawing(Images[(int)ArrayObjects.Header], borders[(int)Object.Header]));
			form.Children.Add(new ImageDrawing(Images[(int)ArrayObjects.Body], borders[(int)Object.Body]));
			form.Children.Add(new ImageDrawing(Images[(int)ArrayObjects.ExitOff], borders[(int)Object.ExitI]));
			form.Children.Add(new ImageDrawing(Images[(int)ArrayObjects.CollapseOff], borders[(int)Object.Collapse]));
			form.Children.Add(new ImageDrawing(Images[(int)ArrayObjects.DropOff], borders[(int)Object.Drop]));
			form.Children.Add(Path);
			form.Children.Add(Rectangle);
			form.Children.Add(TextGD);
			/*form.Children.Add(new ImageDrawing(Images[(int)ArrayObjects.Header], new Rect(0, 0, 500, 30)));
			form.Children.Add(new ImageDrawing(Images[(int)ArrayObjects.Body], new Rect(0, 40, Images[(int)ArrayObjects.Body].Width, Images[(int)ArrayObjects.Body].Height)));
			form.Children.Add(new ImageDrawing(Images[(int)ArrayObjects.ExitOff], new Rect(width - 20 - 30, 5, 20, 20)));
			form.Children.Add(new ImageDrawing(Images[(int)ArrayObjects.CollapseOff], new Rect(width - 60 - 30, 5, 20, 20)));
			form.Children.Add(new ImageDrawing(Images[(int)ArrayObjects.DropOff], new Rect(width - 100 - 30, 5, 20, 20)));*/
		}

		public string Name {
			get { return name; }
			set {
				name = value;
				Update();
			}
		}

		public double Width {
			get { return width; }
			set {
				width = value;
				Update();
			}
		}

		public double Height {
			get { return height; }
			set {
				height = value;
				Update();
			}
		}

		public Rect[] Borders {
			get { return borders; }
			set {
				borders = value;
				Update();
			}
		}

		public Point Coord {
			get { return coord; }
			set {
				coord = value;
				Update();
			}
		}

		public DrawingGroup Form {
			get { return form; }
		}

		bool Test(Rect Border, Point Key) {
			if (Key.X < Border.BottomRight.X && Key.X > Border.TopLeft.X && Key.Y < Border.BottomRight.Y && Key.Y > Border.TopLeft.Y) {
				return true;
			}
			return false;
		}

	}
}
