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
using System.Windows.Shapes;

namespace SoapOpera {
	/// <summary>
	/// Логика взаимодействия для CaptchaForm.xaml
	/// </summary>
	public partial class CaptchaForm : Window {
		public string CaptchaKey;
		public CaptchaForm(Uri Image,  String Name) {
			InitializeComponent();
			this.Title = Name;
			image.Source = new BitmapImage(Image);
		}

		private void button_Click(object sender, RoutedEventArgs e) {
			CaptchaKey = textBox.Text;
			this.Close();
		}
	}
}
