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
using System.Windows.Threading;
using System.Drawing;
using System.Xml;
using System.Runtime.InteropServices;

using System.Windows.Forms;
using System.IO;
using System.Drawing.Imaging;

namespace AutoManager {
	public partial class MainWindow : Window {

		ScriptItem ValidScriptForm() {
			ScriptItem item = new ScriptItem(cbProgram.Text, cbTypeCommand.Text, cbCommand.Text);
			ScriptItem item2 = new ScriptItem(cbProgram.Text, cbTypeCommand.Text, cbCommand.Text);

			while (true) {
				if (cbCounter.IsChecked == true && item2.Command.Contains(options.CounterText) == true) {
					string buff2 = item2.Command.Remove(item2.Command.IndexOf(options.CounterText));
					buff2 = buff2 + 1 + item2.Command.Remove(0, item2.Command.IndexOf(options.CounterText) + options.CounterText.Length);
					item2.Command = buff2;
				}
				else
					break;
			}

			string[] buff = item.Command.ToLower().Split(new Char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
			string[] buff3 = item2.Command.ToLower().Split(new Char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);

			sblInfo.Content = "";

			switch (item.TypeCommand) {
				case "Ожидание":
					try {
						buff3[0] = buff3[0].Replace('.', ',');
						if (Convert.ToDouble(buff3[0]) <= 0) {
							throw new Exception();
						}

						buff[0] = buff[0].Replace(',', '.');
						item.Command = buff[0] + " сек";
					}
					catch {
						sblInfo.Content = "Формат поля Команда \"положительное число сек\" или \"положительное число\" ";
						return null;
					}
					break;
				case "Перемещение мыши":
					try {
						if (buff3[0] == "x:" && buff3[2] == "y:") {
							Convert.ToInt32(buff3[1]);
							Convert.ToInt32(buff3[3]);
							item.Command = "X: " + buff[1] + " Y: " + buff[3];
						}
						else if (buff3[0] == "y:" && buff3[2] == "x:") {
							Convert.ToInt32(buff3[1]);
							Convert.ToInt32(buff3[3]);
							item.Command = "X: " + buff3[3] + " Y: " + buff3[1];
						}
						else {
							Convert.ToInt32(buff3[0]);
							Convert.ToInt32(buff3[1]);
							item.Command = "X: " + buff[0] + " Y: " + buff[1];
						}
					}
					catch {
						sblInfo.Content = "Формат поля Команда \"число число\" или \"x: число y: число\" ";
						return null;
					}
					break;
				case "Повтор":
					try {
						if (buff3.Length > 1 || Convert.ToInt32(buff3[0]) <= 0)
							throw new Exception();
					}
					catch {
						sblInfo.Content = "Формат поля Команда \"число\". Поле не может равнятся 0.";
						return null;
					}
					break;
				default:
					if (item.TypeCommand.Split(' ')[0].Remove(0, item.TypeCommand.IndexOf('.') + 1) == "Edit") {
						item.TypeCommand = item.TypeCommand.Split(' ')[0];
					}
					else if (item.TypeCommand.Split(' ')[0].Remove(0, item.TypeCommand.IndexOf('.') + 1) == "Контекстное") {
						try {
							if (buff3.Length == 3) {
								Convert.ToInt32(buff3[0]);
								Convert.ToInt32(buff3[1]);
								if (Convert.ToInt32(buff3[2]) <= 0)
									throw new Exception();
								item.Command = "X: " + buff[0] + " Y: " + buff[1] + " " + buff[2];
							}
							else if (buff3.Length == 5) {
								if (buff3[0] == "x:" && buff3[2] == "y:") {
									Convert.ToInt32(buff3[1]);
									Convert.ToInt32(buff3[3]);
									if (Convert.ToInt32(buff3[4]) <= 0)
										throw new Exception();
									item.Command = "X: " + buff[1] + " Y: " + buff[3] + " " + buff[4];
								}
								else if (buff3[0] == "y:" && buff3[2] == "x:") {
									Convert.ToInt32(buff3[1]);
									Convert.ToInt32(buff3[3]);
									if (Convert.ToInt32(buff3[4]) <= 0)
										throw new Exception();
									item.Command = "X: " + buff3[3] + " Y: " + buff3[1] + " " + buff[4];
								}
								else
									throw new Exception();
							}
							else
								throw new Exception();
						}
						catch {
							sblInfo.Content = "Формат поля Команда \"число число число\" или \" X: число Y: число число\".";
							return null;
						}
					}
					break;
			}
			return item;
		}

		private IntPtr LowLevelKeyboardHookProc(int nCode, IntPtr wParam, IntPtr lParam) {
			if (nCode < 0) {
				return CppLib.CallNextHookEx(CppLib.m_hHook, nCode, wParam, lParam);
			}
			else {
				var khs = (CppLib.KeyboardHookStruct)Marshal.PtrToStructure(lParam, typeof(CppLib.KeyboardHookStruct));
				if (khs.VirtualKeyCode == (int)VirtualKeyStates.VK_ESCAPE) {
					Stop = true;
				}
				return CppLib.CallNextHookEx(CppLib.m_hHook, nCode, wParam, lParam);
			}
		}

		void UpdatecbProgram() {
			Proc = CppLib.GetProcessList();
			cbProgram.Items.Clear();
			cbCommand.SelectedIndex = -1;
			int OldCount = Proc.Count;
			for (int i = 0; i < OldCount; i++) {
				cbProgram.Items.Add(Proc[i].Name);
				for (int j = 0; j < options.ChildWindows.Count; j++) {
					IntPtr handle;
					if ((handle = CppLib.GetChildWindowHandle(options.ChildWindows[j], Proc[i].Handle)) != IntPtr.Zero) {
						CppLib.Rect rect = new CppLib.Rect();
						CppLib.GetWindowRect(handle, ref rect);
						/*Запомнить на будущее*/
						//Proc.Add(new CppLib.ProcessInfo { Handle = handle, Name = options.ChildWindows[j] + " от " + Proc[i].Name, Border = rect });
						//cbProgram.Items.Add(options.ChildWindows[j] + " от " + Proc[i].Name);
						Proc.Add(new CppLib.ProcessInfo { Handle = handle, Name = options.ChildWindows[j], Border = rect });
						cbProgram.Items.Add(options.ChildWindows[j]);
					}
				}
			}
			cbProgram.Items.Add("Системные команды");
			cbProgram.SelectedIndex = cbProgram.Items.Count - 1;
		}

		void SaveProg() {
			XmlDocument doc = new XmlDocument();
			try {
				doc.Load(Environment.CurrentDirectory + "\\SavePrograms.xml");
			}
			catch {
				doc.AppendChild(doc.CreateElement("Programs"));
			}
			XmlElement progXml = doc.CreateElement("Program");
			progXml.SetAttribute("Name", cbScript.Text);
			doc.ChildNodes[0].AppendChild(progXml);
			for (int i = 0; i < lvScript.Items.Count; i++) {
				XmlElement PoligonXml = doc.CreateElement("Действие_" + i.ToString());
				ScriptItem item = (ScriptItem)lvScript.Items[i];
				PoligonXml.SetAttribute("Программа", item.Program.Replace(' ', '_'));

				PoligonXml.SetAttribute("Тип_команды", item.TypeCommand.Replace(' ', '_'));
				PoligonXml.SetAttribute("Команда", item.Command.Replace(' ', '_'));
				progXml.AppendChild(PoligonXml);
			}
			doc.Save(Environment.CurrentDirectory + "\\SavePrograms.xml");
		}

		void LoadProg(int nom) {
			XmlDocument doc = new XmlDocument();
			try {
				doc.Load(Environment.CurrentDirectory + "\\SavePrograms.xml");
				lvScript.Items.Clear();
				for (int i = 0; i < doc.ChildNodes[0].ChildNodes[nom].ChildNodes.Count; i++) {
					ScriptItem item = new ScriptItem();
					item.Program = doc.ChildNodes[0].ChildNodes[nom].ChildNodes[i].Attributes[0].Value.Replace('_', ' ');
					item.TypeCommand = doc.ChildNodes[0].ChildNodes[nom].ChildNodes[i].Attributes[1].Value.Replace('_', ' ');
					item.Command = doc.ChildNodes[0].ChildNodes[nom].ChildNodes[i].Attributes[2].Value.Replace('_', ' ');
					lvScript.Items.Add(item);
				}
			}
			catch { }
		}

		void DeleteProg(int nom) {
			XmlDocument doc = new XmlDocument();
			try {
				doc.Load(Environment.CurrentDirectory + "\\SavePrograms.xml");
				doc.ChildNodes[0].RemoveChild(doc.ChildNodes[0].ChildNodes[nom]);
				doc.Save(Environment.CurrentDirectory + "\\SavePrograms.xml");
			}
			catch { }
		}

		void LoadcbScript() {
			XmlDocument doc = new XmlDocument();
			try {
				doc.Load(Environment.CurrentDirectory + "\\SavePrograms.xml");
				cbScript.Items.Clear();
				for (int i = 0; i < doc.ChildNodes[0].ChildNodes.Count; i++)
					cbScript.Items.Add(doc.ChildNodes[0].ChildNodes[i].Attributes[0].Value);
				cbScript.Items.Add("Новый сценарий");
				cbScript.Text = "Новый сценарий";
				cbScript.SelectedIndex = cbScript.Items.Count - 1;
			}
			catch { }
		}

		void SaveOption() {
			XmlDocument doc = new XmlDocument();

			doc.AppendChild(doc.CreateElement("Options"));

			XmlElement CounterXml = doc.CreateElement("Counter");
			CounterXml.SetAttribute("Enables", cbCounter.IsChecked.ToString());
			CounterXml.SetAttribute("Text", options.CounterText);
			CounterXml.SetAttribute("Start", options.CounterStart.ToString());
			CounterXml.SetAttribute("Step", options.CounterStep.ToString());
			doc.ChildNodes[0].AppendChild(CounterXml);

			XmlElement ChildwindowsXml = doc.CreateElement("Child_Windows");
			doc.ChildNodes[0].AppendChild(ChildwindowsXml);

			for (int i = 0; i < options.ChildWindows.Count; i++) {
				XmlElement ItemXml = doc.CreateElement("Item");
				ItemXml.SetAttribute("Child_Window", options.ChildWindows[i]);
				ChildwindowsXml.AppendChild(ItemXml);
			}

			XmlElement EllementsXml = doc.CreateElement("Ellements");
			doc.ChildNodes[0].AppendChild(EllementsXml);

			for (int i = 0; i < options.Ellements.Count; i++) {
				XmlElement ItemXml = doc.CreateElement("Item");
				ItemXml.SetAttribute("Ellements", options.Ellements[i]);
				EllementsXml.AppendChild(ItemXml);
			}

			doc.Save(Environment.CurrentDirectory + "\\SaveOptions.xml");
		}

		Options LoadOptions() {
			Options opt = new Options();
			XmlDocument doc = new XmlDocument();
			try {
				doc.Load(Environment.CurrentDirectory + "\\SaveOptions.xml");
				XmlNode CounterXml = doc.ChildNodes[0].ChildNodes[0];

				if (CounterXml.Attributes[0].Value == "True")
					cbCounter.IsChecked = true;
				else
					cbCounter.IsChecked = false;

				try {
					tbCounter.Text = CounterXml.Attributes[1].Value;
				}
				catch { }

				opt.CounterText = CounterXml.Attributes[1].Value;

				try {
					tbCounterStart.Text = CounterXml.Attributes[2].Value;
				}
				catch { }

				opt.CounterStart = Convert.ToInt32(CounterXml.Attributes[2].Value);

				try {
					tbCounterStep.Text = CounterXml.Attributes[3].Value;
				}
				catch { }

				opt.CounterStep = Convert.ToInt32(CounterXml.Attributes[3].Value);

				XmlNode ChildWindowsXml = doc.ChildNodes[0].ChildNodes[1];

				cbChildWindow.Items.Clear();
				for (int i = 0; i < ChildWindowsXml.ChildNodes.Count; i++) {
					opt.ChildWindows.Add(ChildWindowsXml.ChildNodes[i].Attributes[0].Value);
					cbChildWindow.Items.Add(ChildWindowsXml.ChildNodes[i].Attributes[0].Value);
				}

				XmlNode EllementsXml = doc.ChildNodes[0].ChildNodes[2];

				cbEllements.Items.Clear();
				for (int i = 0; i < EllementsXml.ChildNodes.Count; i++) {
					opt.Ellements.Add(EllementsXml.ChildNodes[i].Attributes[0].Value);
					cbEllements.Items.Add(EllementsXml.ChildNodes[i].Attributes[0].Value);
				}

				cbEllements.SelectedIndex = cbChildWindow.SelectedIndex = -1;
			}
			catch { }
			return opt;
		}

		void ValidCounterForm() {
			sblInfo.Content = "";
			if (cbCounter.IsChecked == true) {
				bSaveOption.IsEnabled = true;
				if (tbCounter.Text.Length != 0 && tbCounter.Text.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries).Length == 1) {
					options.CounterText = tbCounter.Text;
				}
				else {
					bSaveOption.IsEnabled = false;
					sblInfo.Content = "Счётчик это одно слово";
				}

				try {
					options.CounterStart = Convert.ToInt32(tbCounterStart.Text);
				}
				catch {
					options.CounterStart = 1;
					sblInfo.Content = "Неправельное начальное значение счётчика. Будет использоваться значение по умолчанию.";
				}

				try {
					options.CounterStep = Convert.ToInt32(tbCounterStep.Text);
				}
				catch {
					options.CounterStep = 1;
					sblInfo.Content = "Неправельный шаг счётчика. Будет использоваться значение по умолчанию.";
				}

			}
			else {
				bSaveOption.IsEnabled = true;
				sblInfo.Content = "";
			}
		}

		void UpdateScript(int j) {
			try {
				XmlDocument doc = new XmlDocument();
				doc.Load(Environment.CurrentDirectory + "\\SavePrograms.xml");
				XmlElement progXml = doc.CreateElement("Program");
				progXml.SetAttribute("Name", cbScript.Text);
				doc.ChildNodes[0].InsertAfter(progXml, doc.ChildNodes[0].ChildNodes[j]);
				doc.ChildNodes[0].RemoveChild(doc.ChildNodes[0].ChildNodes[j]);
				for (int i = 0; i < lvScript.Items.Count; i++) {
					XmlElement PoligonXml = doc.CreateElement("Действие_" + i.ToString());
					ScriptItem item = (ScriptItem)lvScript.Items[i];
					PoligonXml.SetAttribute("Программа", item.Program.Replace(' ', '_'));

					PoligonXml.SetAttribute("Тип_команды", item.TypeCommand.Replace(' ', '_'));
					PoligonXml.SetAttribute("Команда", item.Command.Replace(' ', '_'));
					progXml.AppendChild(PoligonXml);
				}
				doc.Save(Environment.CurrentDirectory + "\\SavePrograms.xml");
			}
			catch { }
		}

		System.Drawing.Color[] GetColorArr(Bitmap picture) {
			System.Drawing.Color[] colors = new System.Drawing.Color[picture.Width * picture.Height];
			byte[] arr;

			using (MemoryStream memoryStream = new MemoryStream()) {
				picture.Save(memoryStream, ImageFormat.Bmp);
				arr = memoryStream.ToArray();
			}

			int w = picture.Height - 1, j = 1;

			for (int i = 54; i < arr.Length; i += 4) {
				System.Drawing.Color color = System.Drawing.Color.FromArgb(arr[i + 3], arr[i + 2], arr[i + 1], arr[i]);
				colors[(picture.Width * w + j) - 1] = color;
				if (j == picture.Width) {
					j = 1;
					w--;
				}
				else
					j++;
			}
			return colors;
		}

		BitmapImage BitmapToImageSource(Bitmap bitmap) {
			using (MemoryStream memory = new MemoryStream()) {
				bitmap.Save(memory, System.Drawing.Imaging.ImageFormat.Bmp);
				memory.Position = 0;
				BitmapImage bitmapimage = new BitmapImage();
				bitmapimage.BeginInit();
				bitmapimage.StreamSource = memory;
				bitmapimage.CacheOption = BitmapCacheOption.OnLoad;
				bitmapimage.EndInit();

				return bitmapimage;
			}
		}
	}
}
