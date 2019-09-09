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
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace prob3 {
	/// <summary>
	/// Логика взаимодействия для MainWindow.xaml
	/// </summary>
	public partial class MainWindow : Window {
		PlayerForm Play;
		Point LeftKey;
		public MainWindow() {
			InitializeComponent();
			Play = new PlayerForm(image,this);			
		}

		private void image_MouseLeftButtonDown(object sender, MouseButtonEventArgs e) {
			LeftKey = e.GetPosition(this);
			Play.Player_MouseDownKey(LeftKey);
		}

		private void image_MouseLeftButtonUp(object sender, MouseButtonEventArgs e) { 
			Play.Player_MouseUpKey(LeftKey);
		}

		private void image_MouseMove(object sender, MouseEventArgs e) {
			LeftKey = e.GetPosition(this);
			Play.Player_MouseMoveKey(LeftKey);
		}

		private void image_KeyDown(object sender, KeyEventArgs e) { }
	}
}
