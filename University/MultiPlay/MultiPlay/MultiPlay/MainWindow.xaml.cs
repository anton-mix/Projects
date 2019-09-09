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
using Microsoft.Win32;
using System.Windows.Forms;
using System.IO;

namespace MultiPlay   {
    public partial class MainWindow : Window  {
        private delegate void UpdateProgressBarDelegate(System.Windows.DependencyProperty dp, Object value);
        List<string> MP3 = new List<string>();
        public MainWindow()  {
            InitializeComponent();
        }
        private void NButton1_Click(object sender, RoutedEventArgs e)  {
            var diolog = new System.Windows.Forms.FolderBrowserDialog();
            if(diolog.ShowDialog() == System.Windows.Forms.DialogResult.OK)
                NTextBox1.Text = diolog.SelectedPath;
        }
        private void RButton1_Click(object sender, RoutedEventArgs e)  {
            var diolog = new System.Windows.Forms.FolderBrowserDialog();
            if (diolog.ShowDialog() == System.Windows.Forms.DialogResult.OK)
                RTextBox1.Text = diolog.SelectedPath;
        }
        private void RButton2_Click(object sender, RoutedEventArgs e)  {
            var diolog = new System.Windows.Forms.FolderBrowserDialog();
            if (diolog.ShowDialog() == System.Windows.Forms.DialogResult.OK)
                RTextBox2.Text = diolog.SelectedPath;
        }
        private void NButton2_Click(object sender, RoutedEventArgs e)  {
            NButton2.IsEnabled = false;
            string[] name = Directory.GetFiles(NTextBox1.Text, "*" + NComboBox.SelectedValue);
			string errorsS = "";
            var random = new Random(DateTime.Now.Millisecond);
			bool errors = false;
            name = name.OrderBy(x => random.Next()).ToArray();
            UpdateProgressBarDelegate updatePbDelegate = new UpdateProgressBarDelegate(NProgressBar1.SetValue);

            for (int i = 0; i < name.Length;i++)  {
                Dispatcher.Invoke(updatePbDelegate, System.Windows.Threading.DispatcherPriority.Background, new object[] {System.Windows.Controls.ProgressBar.ValueProperty, (double)((float)(i + 1) / (float)name.Length * (float)100) });

                string curname = System.IO.Path.GetFileNameWithoutExtension(name[i]);
                string curpath = System.IO.Path.GetDirectoryName(name[i]) + "\\";
                for (int j = 1; j <= 5; j++)
                    if (curname[j] == '.')  {
                        curname = curname.Remove(0, j + 1);
                        break;
                    }
                if(NCheckBox1.IsChecked.Value)
					curname = (i + 1) + "." + curname;
				try {
					var filesmp3 = TagLib.File.Create(name[i]);
					filesmp3.Tag.Title = curname;
					filesmp3.Tag.Track = (uint)(i + 1);
					filesmp3.Tag.TrackCount = (uint)name.Length;
					filesmp3.Save();
					File.Move(name[i],(curpath + curname + NComboBox.SelectedValue));
				}
				catch {
					errorsS = errorsS + curname + "\n";
					errors = true;
				}
            }
            NButton2.IsEnabled = true;
            NProgressBar1.Value = 0;
			if(errors) {
				System.Windows.MessageBox.Show("Errors:\n" + errorsS);
			}
			else {
				System.Windows.MessageBox.Show("Rename complet");
			}
        }
        private void NTextBox1_TextChanged(object sender, TextChangedEventArgs e)  {
            NButton2.IsEnabled = false;
            if (NTextBox1.Text != "")  {
                string[] name = Directory.GetFiles(NTextBox1.Text, "*" + NComboBox.SelectedValue);
                if (name.Length > 1)
                    NButton2.IsEnabled = true;
            }
        }
        private void NComboBox_SelectionChanged(object sender, SelectionChangedEventArgs e)  {
            NButton2.IsEnabled = false;
            if (NTextBox1.Text != "")  {
                string[] name = Directory.GetFiles(NTextBox1.Text, "*" + NComboBox.SelectedValue);
                if (name.Length > 1)
                    NButton2.IsEnabled = true;
            }
        }
        private void RComboBox_SelectionChanged(object sender, SelectionChangedEventArgs e)  {
            RButton3.IsEnabled = false;
            if (RTextBox1.Text != "" && RTextBox2.Text != "")  {
                MP3.Clear();    
                CopyArr(RTextBox1.Text);
                if (MP3.Count > 1)
                    RButton3.IsEnabled = true;
            }
        }
        private void RTextBox1_TextChanged(object sender, TextChangedEventArgs e)  {
            RButton3.IsEnabled = false;
            if (RTextBox1.Text != "" && RTextBox2.Text != "")   {
                MP3.Clear();
                CopyArr(RTextBox1.Text);
                if (MP3.Count > 1)
                    RButton3.IsEnabled = true;
            }
        }
        private void RTextBox2_TextChanged(object sender, TextChangedEventArgs e)  {
            RButton3.IsEnabled = false;
            if (RTextBox1.Text != "" && RTextBox2.Text != "")  {
                MP3.Clear();    
                CopyArr(RTextBox1.Text);
                if (MP3.Count > 1)
                    RButton3.IsEnabled = true;
            }
        }
        private void RButton3_Click(object sender, RoutedEventArgs e)  {
            RButton3.IsEnabled = false;
            MP3.Clear();
            CopyArr(RTextBox1.Text);
            UpdateProgressBarDelegate updatePbDelegate = new UpdateProgressBarDelegate(RProgressBar1.SetValue);
            for (int i = 0; i < MP3.Count; i++)  {
                Dispatcher.Invoke(updatePbDelegate, System.Windows.Threading.DispatcherPriority.Background, new object[] { System.Windows.Controls.ProgressBar.ValueProperty, (double)((float)(i + 1) / (float)MP3.Count * (float)100) });
                File.Copy(MP3[i], RTextBox2.Text + "\\" + System.IO.Path.GetFileName(MP3[i]), true);
            }
            System.Windows.MessageBox.Show("Copy complet");
            RProgressBar1.Value = 0;
            RButton3.IsEnabled = true;
        }
        void CopyArr(string Path)  {
            for (int i = 0; i < Directory.GetDirectories(Path).Length; i++)
                CopyArr(Directory.GetDirectories(Path)[i]);

            for (int i = 0; i < Directory.GetFiles(Path).Length; i++)  {
                string File = Directory.GetFiles(Path)[i];
                bool flag = true;
                for (int j = RComboBox.Text.Length - 1; j >= 0; j--)  {
                    if (RCheckBox1.IsChecked.Value)  {
                        if ((RComboBox.Text).ToLower()[j] != File.ToLower()[File.Length - ((RComboBox.Text).Length - j)])  {
                            flag = false;
                            break;
                        }
                    }
                    else
                        if (((string)RComboBox.SelectedItem)[j] != File[File.Length - (((string)RComboBox.SelectedItem).Length - j)])  {
                            flag = false;
                            break;
                        }
                }
                if (flag)
                    MP3.Add(File);
            }
        }
    }
}
