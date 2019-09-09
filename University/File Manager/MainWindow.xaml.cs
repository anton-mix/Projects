using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Threading;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace prob1 {
	/// <summary>
	/// Логика взаимодействия для MainWindow.xaml
	/// </summary>
	public partial class MainWindow : Window {
		MyStruct prob;
		Node prob2, prob3, prob4;
		TranslateTransform Trans = new TranslateTransform();
		DispatcherTimer timer = new DispatcherTimer();
		public MainWindow() {
			InitializeComponent();
			prob = new MyStruct(MyCanvas,0,0);
			prob.TransCenterX = MyCanvas.Width / 2;
			prob.TransCenterY = MyCanvas.Height / 2;
			timer.Tick += new EventHandler(timer_Tick);
			timer.Interval = TimeSpan.FromSeconds(5);
			//prob.InitFerstStruct();
			//timer.Start();
			prob2 = new Node(-(int)MyCanvas.Width / 2,-(int)MyCanvas.Height / 2,15,1,MyCanvas,1,Trans);
			Trans.X = MyCanvas.Width;
			Trans.Y = MyCanvas.Height;
			//prob2.DrawLine((LineDir)0,75,0,0);
			prob2.DrawLine((LineDir)1,75,0,0);
			//prob2.DrawLine((LineDir)2,75,0,0);
			prob2.DrawLine((LineDir)3,75,0,0);
			//prob2.DrawLine((LineDir)4,75,0,0);
			prob2.DrawLine((LineDir)5,75,0,0);
			//prob2.DrawLine((LineDir)6,75,0,0);
			prob2.DrawLine((LineDir)7,75,0,0);
			//prob2.MoveNode(50,50);*/
		}

		private void BackGround_MouseDown(object sender,MouseButtonEventArgs e) {
			this.DragMove();
		}

		private void button_Click(object sender,RoutedEventArgs e) {
			prob.deleteStruct();
			prob.InitFerstStruct();
		}

		private void button1_Click(object sender,RoutedEventArgs e) {
			this.Close();
		}

		private void button2_Click(object sender,RoutedEventArgs e) {
			prob.PassivAnimation();
		}

		private void textBox_KeyUp(object sender,KeyEventArgs e) {
			double duff;
			if (Double.TryParse(textBox.Text,out duff))
				prob3.Center = new Point(duff,prob3.Center.Y);
		}

		private void textBox1_KeyUp(object sender,KeyEventArgs e) {
			double duff;
			if (Double.TryParse(textBox1.Text,out duff))
				prob3.Center = new Point(prob3.Center.X,duff);
		}
		private void timer_Tick(object sender,EventArgs e) {
			prob.PassivAnimation();
		}
	}
}
